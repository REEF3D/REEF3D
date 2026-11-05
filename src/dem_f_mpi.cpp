/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
Domain decomposition of the DEM.

Tier 0 (distributed): a particle is owned by the rank whose subdomain holds its
centroid. Copies (ghosts) are sent to every neighbouring rank whose subdomain is
within the particle's halo: contact range plus the reach of its fluid coupling.
Ownership changes after each DEM substep (migration).

Tier 1 (replicated): particles too large for a halo within one subdomain (E 24)
are kept on all ranks. Their contact impulses and fluid contributions are summed
with global reductions over the (few) tier 1 particles only.

Contacts are solved by the rank that owns the lower global id (tier 1 - tier 1:
rank 0; tier 0 - tier 1: owner of the tier 0 particle). Contacts of particles
solved on one rank only are swept locally; contacts of particles shared by
several ranks are swept in colour phases (one colour at a time), each followed
by returning the ghost velocity corrections to the owners and sending the new
velocities back, so the iteration stays Gauss-Seidel across rank boundaries.
--------------------------------------------------------------------*/

#include"dem_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<algorithm>
#include<cstring>

// ---------------------------------------------------------------------------------------------
// setup: subdomain boxes, tiers, partner ranks, initial distribution
// ---------------------------------------------------------------------------------------------

void dem_f::decomp_ini(lexer *p, ghostcell *pgc)
{
    myrank = p->mpirank;
    nproc = p->mpi_size;
    core.myrank = myrank;

    // subdomain boxes of all ranks (inactive directions unbounded)
    double mybox[6] = {p->originx,p->endx,p->originy,p->endy,p->originz,p->endz};
    if(p->j_dir==0)
    {
        mybox[2] = -1.0e30;
        mybox[3] =  1.0e30;
    }
    if(solver==5)
    {
        mybox[4] = -1.0e30;
        mybox[5] =  1.0e30;
    }
    p_gmax[0] = p->global_xmax;
    p_gmax[1] = p->j_dir==1 ? p->global_ymax : 1.0e30;
    p_gmax[2] = solver==6 ? p->global_zmax : 1.0e30;

    boxes.assign(6*nproc,0.0);
    MPI_Allgather(mybox,6,MPI_DOUBLE,boxes.data(),6,MPI_DOUBLE,pgc->mpi_comm);

    // smallest subdomain extent in the active directions
    double minext = 1.0e30;
    for(int r=0; r<nproc; ++r)
    for(int d=0; d<3; ++d)
    {
        double ext = boxes[6*r+2*d+1] - boxes[6*r+2*d];
        if(ext<1.0e20)
        minext = std::min(minext,ext);
    }

    // tiers: distributed if the particle is small compared with the subdomains
    int n0=0, n1=0;
    rmax0 = 0.0;
    for(auto &B : core.bodies)
    {
        double rb = core.shapes[B.shape].rbound;
        B.ghost = false;

        if(nproc==1 || (p->E24>0.0 && rb<=p->E24*minext))
        {
            B.tier = 0;
            rmax0 = std::max(rmax0,rb);
            ++n0;
        }
        else
        {
            B.tier = 1;
            B.owner = -1;
            ++n1;
        }
    }

    // halo reach and partner ranks: all ranks whose subdomain is closer than the largest halo
    double hmax = 0.0;
    for(auto &B : core.bodies)
    if(B.tier==0)
    hmax = std::max(hmax,halo(B));
    hmax += dxs;

    partners.clear();
    for(int r=0; r<nproc; ++r)
    {
        if(r==myrank)
        continue;

        double d2=0.0;
        for(int d=0; d<3; ++d)
        {
            double lo = std::max(mybox[2*d],boxes[6*r+2*d]);
            double hi = std::min(mybox[2*d+1],boxes[6*r+2*d+1]);
            if(lo>hi)
            d2 += (lo-hi)*(lo-hi);
        }
        if(sqrt(d2)<hmax)
        partners.push_back(r);
    }

    partner_index.assign(nproc,-1);
    for(size_t k=0; k<partners.size(); ++k)
    partner_index[partners[k]] = k;

    // colours for the shared-contact phases: ranks closer than twice the halo may share a particle
    // and get different colours (greedy in rank order, identical on all ranks)
    {
        vector<int> col(nproc,-1);
        ncolors = 1;
        for(int r=0; r<nproc; ++r)
        {
            vector<bool> used(nproc+1,false);
            for(int s=0; s<r; ++s)
            {
                double d2=0.0;
                for(int d=0; d<3; ++d)
                {
                    double lo = std::max(boxes[6*r+2*d],boxes[6*s+2*d]);
                    double hi = std::min(boxes[6*r+2*d+1],boxes[6*s+2*d+1]);
                    if(lo>hi)
                    d2 += (lo-hi)*(lo-hi);
                }
                if(sqrt(d2)<2.0*hmax)
                used[col[s]] = true;
            }
            int c=0;
            while(used[c])
            ++c;
            col[r] = c;
            ncolors = std::max(ncolors,c+1);
        }
        mycolor = col[myrank];
    }

    // keep the replicated particles and the owned distributed ones; order: tier 1 first
    vector<dem_body> keep;
    int lost=0;
    for(auto &B : core.bodies)
    if(B.tier==1)
    keep.push_back(B);

    for(auto &B : core.bodies)
    if(B.tier==0)
    {
        int r = find_owner(B.x);
        if(r<0)
        {
            ++lost;
            continue;
        }
        if(r==myrank)
        {
            B.owner = myrank;
            keep.push_back(B);
        }
    }
    core.bodies.swap(keep);
    nb = core.bodies.size();

    if(myrank==0)
    {
        cout<<"DEM: domain decomposition over "<<nproc<<" ranks: "<<n0<<" distributed particles, "<<n1<<" replicated particles";
        if(nproc>1)
        cout<<", halo "<<hmax<<" m, colours "<<ncolors;
        cout<<endl;
        if(lost>0)
        cout<<"DEM: warning, "<<lost<<" particles are outside all subdomains and are removed"<<endl;
    }

    int npmax = partners.size();
    MPI_Allreduce(MPI_IN_PLACE,&npmax,1,MPI_INT,MPI_MAX,pgc->mpi_comm);
    if(myrank==0 && nproc>1)
    cout<<"DEM: partner ranks per subdomain (max): "<<npmax<<endl;

    ghost_exchange(p,pgc);
}

double dem_f::halo(const dem_body &B)
{
    // contact range, resolved forcing and internal momentum, unresolved kernel, NHFLOW ring
    const dem_shape &S = core.shapes[B.shape];
    double rb = S.rbound;
    double h = rb + rmax0 + std::max(rmax0,std::max(2.0*dxs,2.0*travel*rmax0));
    h = std::max(h, std::max(std::max(S.deq,void_factor*S.deq),kernel_cells*dxs) + 2.0*dxs);
    h = std::max(h, rb + (hs_factor+3.0)*dxs);
    return h;
}

bool dem_f::in_box(int r, const dem_vec &x)
{
    // consistent with owns(): half-open intervals, upper global boundary included
    auto in = [](double v, double lo, double hi, double gmax)
    {
        return v>=lo && (v<hi || (hi>=gmax-1.0e-12 && v<=hi));
    };
    const double *b = &boxes[6*r];
    return in(x(0),b[0],b[1],p_gmax[0]) && in(x(1),b[2],b[3],p_gmax[1]) && in(x(2),b[4],b[5],p_gmax[2]);
}

int dem_f::find_owner(const dem_vec &x)
{
    if(nproc==1)
    return 0;

    if(in_box(myrank,x))
    return myrank;

    for(int r : partners)
    if(in_box(r,x))
    return r;

    for(int r=0; r<nproc; ++r)
    if(in_box(r,x))
    return r;

    return -1;
}

double dem_f::boxdist(int r, const dem_vec &x)
{
    const double *b = &boxes[6*r];
    double d2=0.0;
    for(int d=0; d<3; ++d)
    {
        double e = 0.0;
        if(x(d)<b[2*d]) e = b[2*d]-x(d);
        else if(x(d)>b[2*d+1]) e = x(d)-b[2*d+1];
        d2 += e*e;
    }
    return sqrt(d2);
}

// ---------------------------------------------------------------------------------------------
// point-to-point exchange with the partner ranks
// ---------------------------------------------------------------------------------------------

void dem_f::p2p(ghostcell *pgc, vector<vector<double>> &send, vector<vector<double>> &recv)
{
    int np = partners.size();
    recv.assign(np,vector<double>());
    if(np==0)
    return;

    vector<long long> scount(np), rcount(np);
    vector<MPI_Request> req(2*np);

    for(int k=0; k<np; ++k)
    {
        scount[k] = send[k].size();
        MPI_Irecv(&rcount[k],1,MPI_LONG_LONG,partners[k],7301,pgc->mpi_comm,&req[k]);
        MPI_Isend(&scount[k],1,MPI_LONG_LONG,partners[k],7301,pgc->mpi_comm,&req[np+k]);
    }
    MPI_Waitall(2*np,req.data(),MPI_STATUSES_IGNORE);

    int nreq=0;
    for(int k=0; k<np; ++k)
    {
        recv[k].resize(rcount[k]);
        if(rcount[k]>0)
        MPI_Irecv(recv[k].data(),int(rcount[k]),MPI_DOUBLE,partners[k],7302,pgc->mpi_comm,&req[nreq++]);
        if(scount[k]>0)
        MPI_Isend(send[k].data(),int(scount[k]),MPI_DOUBLE,partners[k],7302,pgc->mpi_comm,&req[nreq++]);
    }
    MPI_Waitall(nreq,req.data(),MPI_STATUSES_IGNORE);
}

void dem_f::rebuild_maps()
{
    nb = core.bodies.size();
    gid2loc.clear();
    tier1idx.clear();
    for(int n=0; n<nb; ++n)
    {
        gid2loc[core.bodies[n].id] = n;
        if(core.bodies[n].tier==1)
        tier1idx.push_back(n);
    }
    std::sort(tier1idx.begin(),tier1idx.end(),[&](int a, int b){return core.bodies[a].id<core.bodies[b].id;});
}

// ---------------------------------------------------------------------------------------------
// ghosts
// ---------------------------------------------------------------------------------------------

void dem_f::erase_ghosts()
{
    size_t k=0;
    for(size_t n=0; n<core.bodies.size(); ++n)
    if(!core.bodies[n].ghost)
    {
        if(k!=n)
        core.bodies[k] = core.bodies[n];
        ++k;
    }
    core.bodies.resize(k);
    ghostdest.clear();
    rebuild_maps();
}

void dem_f::ghost_exchange(lexer *p, ghostcell *pgc)
{
    erase_ghosts();

    int np = partners.size();
    ghostdest.assign(core.bodies.size(),vector<int>());

    if(np>0)
    {
        vector<vector<double>> send(np), recv;

        for(size_t n=0; n<core.bodies.size(); ++n)
        {
            const dem_body &B = core.bodies[n];
            if(B.tier!=0 || !B.active)
            continue;

            double h = halo(B);
            for(int k=0; k<np; ++k)
            if(boxdist(partners[k],B.x)<h)
            {
                ghostdest[n].push_back(k);
                dem_writer w(send[k]);
                w.d(B.id); w.d(B.shape); w.d(B.mat); w.d(B.owner); w.d(B.mode); w.d(B.cpl.basemode);
                w.d(B.fixed?1:0); w.d(B.m); w.v(B.Ib);
                w.v(B.x); w.q(B.q); w.v(B.v); w.v(B.w); w.v(B.vold); w.v(B.wold);
                w.v(B.invm); w.m(B.invI);
            }
        }

        p2p(pgc,send,recv);

        for(int k=0; k<np; ++k)
        {
            dem_reader r(recv[k],0);
            while(r.pos<recv[k].size())
            {
                dem_body B;
                B.id = r.i(); B.shape = r.i(); B.mat = r.i(); B.owner = r.i(); B.mode = r.i(); B.cpl.basemode = r.i();
                B.fixed = r.i()==1; B.m = r.d(); B.Ib = r.v();
                B.x = r.v(); B.q = r.q(); B.v = r.v(); B.w = r.v(); B.vold = r.v(); B.wold = r.v();
                B.invm = r.v(); B.invI = r.m();
                B.R = B.q.toRotationMatrix();
                B.tier = 0;
                B.ghost = true;
                B.active = true;
                core.bodies.push_back(B);
            }
        }
    }

    rebuild_maps();
    ghostdest.resize(core.bodies.size());

    // reference velocities for the solver synchronisation
    for(auto &B : core.bodies)
    {
        B.vs = B.v;
        B.ws = B.w;
    }
}

// ---------------------------------------------------------------------------------------------
// migration of owned particles to the rank that holds their centroid
// ---------------------------------------------------------------------------------------------

void dem_f::pack_body(vector<double> &buf, const dem_body &B)
{
    dem_writer w(buf);
    const dem_cpl &c = B.cpl;
    w.d(B.id); w.d(B.shape); w.d(B.mat); w.d(B.mode); w.d(B.fixed?1:0); w.d(B.active?1:0);
    w.d(B.m); w.v(B.Ib); w.v(B.x); w.q(B.q); w.v(B.v); w.v(B.w); w.v(B.vold); w.v(B.wold);
    w.v(B.F); w.v(B.T); w.d(B.K); w.d(B.Kr); w.d(B.madd); w.v(B.uf); w.v(B.af);

    w.d(c.rhof); w.d(c.nuf); w.d(c.epsf); w.d(c.vsub); w.v(c.ufl); w.v(c.ufl_old); w.v(c.Fb); w.v(c.Tb);
    w.d(c.fluidcount); w.d(c.ufl_valid?1:0); w.v(c.Fibm); w.v(c.Tibm); w.d(c.mfl); w.d(c.hvol);
    for(int s=0; s<3; ++s) {w.v(c.Fs[s]); w.v(c.Ts[s]);}
    w.v(c.Ifl); w.v(c.Lfl); w.v(c.Ifl_old); w.v(c.Lfl_old); w.d(c.Ifl_valid?1:0);
    w.v(c.vprev); w.v(c.wprev); w.v(c.Ffp); w.v(c.Fhyd);
    for(int s=0; s<4; ++s) w.d(c.sw[s]);
    w.d(c.basemode);
}

void dem_f::unpack_body(dem_reader &r, dem_body &B)
{
    dem_cpl &c = B.cpl;
    B.id = r.i(); B.shape = r.i(); B.mat = r.i(); B.mode = r.i(); B.fixed = r.i()==1; B.active = r.i()==1;
    B.m = r.d(); B.Ib = r.v(); B.x = r.v(); B.q = r.q(); B.v = r.v(); B.w = r.v(); B.vold = r.v(); B.wold = r.v();
    B.F = r.v(); B.T = r.v(); B.K = r.d(); B.Kr = r.d(); B.madd = r.d(); B.uf = r.v(); B.af = r.v();

    c.rhof = r.d(); c.nuf = r.d(); c.epsf = r.d(); c.vsub = r.d(); c.ufl = r.v(); c.ufl_old = r.v(); c.Fb = r.v(); c.Tb = r.v();
    c.fluidcount = r.i(); c.ufl_valid = r.i()==1; c.Fibm = r.v(); c.Tibm = r.v(); c.mfl = r.d(); c.hvol = r.d();
    for(int s=0; s<3; ++s) {c.Fs[s] = r.v(); c.Ts[s] = r.v();}
    c.Ifl = r.v(); c.Lfl = r.v(); c.Ifl_old = r.v(); c.Lfl_old = r.v(); c.Ifl_valid = r.i()==1;
    c.vprev = r.v(); c.wprev = r.v(); c.Ffp = r.v(); c.Fhyd = r.v();
    for(int s=0; s<4; ++s) c.sw[s] = r.d();
    c.basemode = r.i();

    B.R = B.q.toRotationMatrix();
    B.tier = 0;
    B.ghost = false;
}

void dem_f::migrate(lexer *p, ghostcell *pgc)
{
    erase_ghosts();

    if(nproc==1)
    return;

    int np = partners.size();
    vector<vector<double>> send(np), recv;
    vector<dem_body> keep;
    keep.reserve(core.bodies.size());
    int far=0;

    for(auto &B : core.bodies)
    {
        if(B.tier==1 || !B.active)
        {
            keep.push_back(B);
            continue;
        }

        int r = find_owner(B.x);
        if(r<0 || r==myrank)
        {
            keep.push_back(B);      // outside all subdomains: deactivate() handles it
            continue;
        }

        int k = partner_index[r];
        if(k<0)
        {
            ++far;
            keep.push_back(B);
            continue;
        }
        pack_body(send[k],B);
    }

    core.bodies.swap(keep);
    p2p(pgc,send,recv);

    for(int k=0; k<np; ++k)
    {
        dem_reader r(recv[k],0);
        while(r.pos<recv[k].size())
        {
            dem_body B;
            unpack_body(r,B);
            B.owner = myrank;
            B.vs = B.v;
            B.ws = B.w;
            core.bodies.push_back(B);
        }
    }

    if(far>0)
    cout<<"DEM: warning, rank "<<myrank<<": "<<far<<" particles moved beyond the partner ranks and stay with their owner"<<endl;

    rebuild_maps();
}

// ---------------------------------------------------------------------------------------------
// reductions: ghost contributions to the owner, replicated particles over all ranks
// ---------------------------------------------------------------------------------------------

void dem_f::reduce_owner(ghostcell *pgc, vector<double> &buf, int nv, bool toghosts, bool usemin)
{
    int np = partners.size();

    // ghosts -> owners
    if(np>0)
    {
        vector<vector<double>> send(np), recv;
        for(int n=0; n<nb; ++n)
        {
            const dem_body &B = core.bodies[n];
            if(!B.ghost)
            continue;
            int k = partner_index[B.owner];
            if(k<0)
            continue;
            send[k].push_back(B.id);
            for(int q=0; q<nv; ++q)
            send[k].push_back(buf[n*nv+q]);
        }

        p2p(pgc,send,recv);

        for(int k=0; k<np; ++k)
        for(size_t pos=0; pos+nv+1<=recv[k].size(); pos+=nv+1)
        {
            auto it = gid2loc.find(int(llround(recv[k][pos])));
            if(it==gid2loc.end())
            continue;
            int n = it->second;
            for(int q=0; q<nv; ++q)
            {
                double v = recv[k][pos+1+q];
                buf[n*nv+q] = usemin ? std::min(buf[n*nv+q],v) : buf[n*nv+q]+v;
            }
        }
    }

    // replicated particles
    int n1 = tier1idx.size();
    if(n1>0 && nproc>1)
    {
        vector<double> t(n1*nv);
        for(int m=0; m<n1; ++m)
        for(int q=0; q<nv; ++q)
        t[m*nv+q] = buf[tier1idx[m]*nv+q];

        MPI_Allreduce(MPI_IN_PLACE,t.data(),n1*nv,MPI_DOUBLE,usemin ? MPI_MIN : MPI_SUM,pgc->mpi_comm);

        for(int m=0; m<n1; ++m)
        for(int q=0; q<nv; ++q)
        buf[tier1idx[m]*nv+q] = t[m*nv+q];
    }

    if(toghosts)
    owner_to_ghosts(pgc,buf,nv);
}

void dem_f::owner_to_ghosts(ghostcell *pgc, vector<double> &buf, int nv)
{
    int np = partners.size();
    if(np==0)
    return;

    vector<vector<double>> send(np), recv;
    for(int n=0; n<nb; ++n)
    {
        if(core.bodies[n].ghost || n>=int(ghostdest.size()))
        continue;
        for(int k : ghostdest[n])
        {
            send[k].push_back(core.bodies[n].id);
            for(int q=0; q<nv; ++q)
            send[k].push_back(buf[n*nv+q]);
        }
    }

    p2p(pgc,send,recv);

    for(int k=0; k<np; ++k)
    for(size_t pos=0; pos+nv+1<=recv[k].size(); pos+=nv+1)
    {
        auto it = gid2loc.find(int(llround(recv[k][pos])));
        if(it==gid2loc.end() || !core.bodies[it->second].ghost)
        continue;
        int n = it->second;
        for(int q=0; q<nv; ++q)
        buf[n*nv+q] = recv[k][pos+1+q];
    }
}

// ---------------------------------------------------------------------------------------------
// contact solver synchronisation
// ---------------------------------------------------------------------------------------------

double dem_f::sync_solver(ghostcell *pgc, int mode, double res, bool final)
{
    // velocity corrections accumulated since the last synchronisation
    vector<double> dv(6*nb,0.0);
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        dem_vec a = mode==0 ? dem_vec(B.v - B.vs) : dem_vec(B.vp - B.vps);
        dem_vec b = mode==0 ? dem_vec(B.w - B.ws) : dem_vec(B.wp - B.wps);
        for(int q=0; q<3; ++q)
        {
            dv[6*n+q] = a(q);
            dv[6*n+3+q] = b(q);
        }
    }

    // sum at the owners (ghost corrections) and over all ranks (replicated particles)
    reduce_owner(pgc,dv,6,false,false);

    // owners: reference plus the corrections of all copies (within a colour phase only one rank
    // modifies a shared distributed particle, so the sum is exact); replicated particles: mean of the
    // corrections (mass splitting, see mass_split)
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(B.ghost)
        continue;
        double f = B.tier==1 ? 1.0/std::max(1.0,B.split) : 1.0;
        dem_vec a = f*dem_vec(dv[6*n],dv[6*n+1],dv[6*n+2]);
        dem_vec b = f*dem_vec(dv[6*n+3],dv[6*n+4],dv[6*n+5]);
        if(mode==0)
        {
            B.v = B.vs + a;
            B.w = B.ws + b;
        }
        else
        {
            B.vp = B.vps + a;
            B.wp = B.wps + b;
        }
    }

    // owners send the corrected velocities to their ghosts
    vector<double> vel(6*nb,0.0);
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        dem_vec a = mode==0 ? B.v : B.vp;
        dem_vec b = mode==0 ? B.w : B.wp;
        for(int q=0; q<3; ++q)
        {
            vel[6*n+q] = a(q);
            vel[6*n+3+q] = b(q);
        }
    }
    owner_to_ghosts(pgc,vel,6);

    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        dem_vec a(vel[6*n],vel[6*n+1],vel[6*n+2]);
        dem_vec b(vel[6*n+3],vel[6*n+4],vel[6*n+5]);
        if(mode==0)
        {
            if(B.ghost)
            {
                B.v = a;
                B.w = b;
            }
            B.vs = B.v;
            B.ws = B.w;
        }
        else
        {
            if(B.ghost)
            {
                B.vp = a;
                B.wp = b;
            }
            B.vps = B.vp;
            B.wps = B.wp;
        }
    }

    if(final)
    MPI_Allreduce(MPI_IN_PLACE,&res,1,MPI_DOUBLE,MPI_MAX,pgc->mpi_comm);
    return res;
}

// ---------------------------------------------------------------------------------------------
// sharing count: number of ranks that solve contacts of a particle; contacts of shared particles
// are swept in the colour phases; replicated particles are mass-split
// ---------------------------------------------------------------------------------------------

void dem_f::mass_split(ghostcell *pgc)
{
    vector<double> f(nb,0.0);
    for(auto &C : core.contacts)
    {
        f[C.a] = 1.0;
        if(C.b>=0)
        f[C.b] = 1.0;
    }

    reduce_owner(pgc,f,1,true,false);

    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        B.split = std::max(1.0,f[n]);

        // replicated particles can be corrected by two ranks of the same colour in one phase:
        // mass splitting (Tonge et al. 2012), each copy is solved with mass m/split and the
        // synchronisation averages the corrections, which keeps the velocity equal to v + invm*sum(P)
        if(B.tier==1 && B.split>1.5)
        {
            B.invm *= B.split;
            B.invI *= B.split;
        }
    }
}

// ---------------------------------------------------------------------------------------------
// wall contacts: routed to the rank that solves them
// ---------------------------------------------------------------------------------------------

void dem_f::route_walls(lexer *p, ghostcell *pgc, const vector<double> &loc, vector<dem_contact> &cts)
{
    // record: local body index, feature, x(3), n(3), gap
    const int nv = 9;
    int np = partners.size();
    vector<vector<double>> send(np), recv;
    vector<double> tier1rec;

    auto make = [&](int n, const double *b)
    {
        dem_contact C;
        C.a = n;
        C.b = -1;
        C.key = core.make_key(core.bodies[n].id,-1,int(llround(b[0])));
        C.x = dem_vec(b[1],b[2],b[3]);
        C.n = dem_vec(b[4],b[5],b[6]);
        C.gap = b[7];
        const dem_material &M = core.mats[core.bodies[n].mat];
        C.mu = 0.5*(M.friction + core.wallmat.friction);
        C.e  = 0.5*(M.restitution + core.wallmat.restitution);
        cts.push_back(C);
    };

    for(size_t c=0; c<loc.size()/nv; ++c)
    {
        const double *b = &loc[c*nv];
        int n = int(llround(b[0]));
        const dem_body &B = core.bodies[n];

        if(B.tier==1)
        {
            tier1rec.push_back(B.id);
            tier1rec.insert(tier1rec.end(),b+1,b+nv);
        }
        else if(!B.ghost)
        make(n,b+1);
        else
        {
            int k = partner_index[B.owner];
            if(k<0)
            continue;
            send[k].push_back(B.id);
            send[k].insert(send[k].end(),b+1,b+nv);
        }
    }

    // ghost records to the owners
    p2p(pgc,send,recv);
    for(int k=0; k<np; ++k)
    for(size_t pos=0; pos<recv[k].size(); pos+=nv)
    {
        auto it = gid2loc.find(int(llround(recv[k][pos])));
        if(it!=gid2loc.end() && !core.bodies[it->second].ghost)
        make(it->second,&recv[k][pos+1]);
    }

    // replicated particles: rank 0 solves their wall contacts
    if(nproc==1)
    {
        for(size_t pos=0; pos<tier1rec.size(); pos+=nv)
        {
            auto it = gid2loc.find(int(llround(tier1rec[pos])));
            if(it!=gid2loc.end())
            make(it->second,&tier1rec[pos+1]);
        }
        return;
    }

    int nloc = tier1rec.size();
    vector<int> counts(nproc,0), displs(nproc,0);
    MPI_Gather(&nloc,1,MPI_INT,counts.data(),1,MPI_INT,0,pgc->mpi_comm);

    int ntot=0;
    if(myrank==0)
    for(int r=0; r<nproc; ++r)
    {
        displs[r] = ntot;
        ntot += counts[r];
    }

    vector<double> all(std::max(ntot,1));
    MPI_Gatherv(tier1rec.data(),nloc,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,0,pgc->mpi_comm);

    if(myrank==0)
    for(int pos=0; pos<ntot; pos+=nv)
    {
        auto it = gid2loc.find(int(llround(all[pos])));
        if(it!=gid2loc.end())
        make(it->second,&all[pos+1]);
    }
}

// ---------------------------------------------------------------------------------------------
// replicated particles: rank 0 broadcasts their state (removes round-off drift)
// ---------------------------------------------------------------------------------------------

void dem_f::sync(lexer *p, ghostcell *pgc)
{
    int n1 = tier1idx.size();
    if(n1==0 || nproc==1)
    return;

    const int nv = 14;
    vector<double> buf(n1*nv);

    if(myrank==0)
    for(int m=0; m<n1; ++m)
    {
        const dem_body &B = core.bodies[tier1idx[m]];
        double *b = &buf[m*nv];
        b[0]=B.x(0); b[1]=B.x(1); b[2]=B.x(2);
        b[3]=B.q.w(); b[4]=B.q.x(); b[5]=B.q.y(); b[6]=B.q.z();
        b[7]=B.v(0); b[8]=B.v(1); b[9]=B.v(2);
        b[10]=B.w(0); b[11]=B.w(1); b[12]=B.w(2);
        b[13]=B.active ? 1.0 : 0.0;
    }

    MPI_Bcast(buf.data(),n1*nv,MPI_DOUBLE,0,pgc->mpi_comm);

    if(myrank>0)
    for(int m=0; m<n1; ++m)
    {
        dem_body &B = core.bodies[tier1idx[m]];
        const double *b = &buf[m*nv];
        B.x = dem_vec(b[0],b[1],b[2]);
        B.q = dem_quat(b[3],b[4],b[5],b[6]);
        B.v = dem_vec(b[7],b[8],b[9]);
        B.w = dem_vec(b[10],b[11],b[12]);
        B.active = b[13]>0.5;
        B.R = B.q.toRotationMatrix();
    }
}

// ---------------------------------------------------------------------------------------------
// global quantities
// ---------------------------------------------------------------------------------------------

void dem_f::global_stats(ghostcell *pgc, double &vmax, double &rmin)
{
    vmax = 0.0;
    rmin = 1.0e20;
    for(auto &B : core.bodies)
    {
        // travel limit per substep: moving particles only, fixed ones are caught by the speculative contacts
        if(B.ghost || !B.active || B.fixed)
        continue;
        rmin = std::min(rmin,core.shapes[B.shape].rbound);
        vmax = std::max(vmax, B.v.norm() + core.shapes[B.shape].rbound*B.w.norm());
    }
    double buf[2] = {vmax,-rmin};
    MPI_Allreduce(MPI_IN_PLACE,buf,2,MPI_DOUBLE,MPI_MAX,pgc->mpi_comm);
    vmax = buf[0];
    rmin = -buf[1];
}

// gather the owned particles on rank 0 for output (sorted by global id)
void dem_f::gather_output(ghostcell *pgc, vector<dem_body> &out)
{
    const int nv = 23;
    vector<double> loc;
    for(auto &B : core.bodies)
    {
        if(B.ghost || (B.tier==1 && myrank!=0))
        continue;
        loc.push_back(B.id); loc.push_back(B.shape); loc.push_back(B.mode); loc.push_back(B.active?1:0);
        for(int q=0; q<3; ++q) loc.push_back(B.x(q));
        loc.push_back(B.q.w()); loc.push_back(B.q.x()); loc.push_back(B.q.y()); loc.push_back(B.q.z());
        for(int q=0; q<3; ++q) loc.push_back(B.v(q));
        for(int q=0; q<3; ++q) loc.push_back(B.w(q));
        for(int q=0; q<3; ++q) loc.push_back(B.cpl.Fhyd(q));
        for(int q=0; q<3; ++q) loc.push_back(B.cpl.ufl(q));
    }

    int nloc = loc.size();
    vector<int> counts(nproc,0), displs(nproc,0);
    MPI_Gather(&nloc,1,MPI_INT,counts.data(),1,MPI_INT,0,pgc->mpi_comm);
    int ntot=0;
    if(myrank==0)
    for(int r=0; r<nproc; ++r)
    {
        displs[r] = ntot;
        ntot += counts[r];
    }
    vector<double> all(std::max(ntot,1));
    MPI_Gatherv(loc.data(),nloc,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,0,pgc->mpi_comm);

    out.clear();
    if(myrank!=0)
    return;

    for(int pos=0; pos+nv<=ntot; pos+=nv)
    {
        dem_reader r(all,pos);
        dem_body B;
        B.id = r.i(); B.shape = r.i(); B.mode = r.i(); B.active = r.i()==1;
        B.x = r.v(); B.q = r.q(); B.v = r.v(); B.w = r.v();
        B.cpl.Fhyd = r.v(); B.cpl.ufl = r.v();
        B.R = B.q.toRotationMatrix();
        out.push_back(B);
    }
    std::sort(out.begin(),out.end(),[](const dem_body &a, const dem_body &b){return a.id<b.id;});
}
