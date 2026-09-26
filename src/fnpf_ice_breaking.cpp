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

#include"fnpf_ice.h"
#include"lexer.h"
#include"ghostcell.h"
#include"ice_contact.h"
#include<iomanip>
#include<algorithm>

// Ice breaking (A 390): a floe splits in two along a straight cut, at most once per floe and check.
//
// Flexural failure (A 390 1,3): internal bending moment of the rigid floe about candidate cut lines
// (A 395: directions x offsets), from the loads on one side of the cut:
//     q = phi*A*( p_lid - rho_i*g*h - rho_i*h*a_z(x,y) ),   a_z = a_G,z + (alpha x r)_z   (d'Alembert)
//     M = sum_{n.r > s} q*(n.r - s),     sigma = 6|M|/(h^2 L),   L = chord of the cut
// q sums to zero over a floe at rest and in rigid-body equilibrium, so there is no spurious moment.
// Breaks when sigma >= sigma_f (A 391) on the cut with the largest stress.
// Rigid-floe statics: an upper bound for floes longer than the flexural length.
//
// Contact splitting (A 390 2,3, 3D): a floe splits along the contact normal through the contact point
// when the pair normal force exceeds F_s = C*K_IC*h*sqrt(D) (A 392), D = 2*sqrt(A/pi); LEFM scaling of
// in-plane splitting (Bhat 1988, Lu et al. 2015), C to be calibrated. The force is the pair normal
// force of the non-smooth contact averaged over t_c (A 396), see contact_average: impact impulses are
// resolved within one or two steps, impulse/dt alone depends on dt and on the timing of the impact.
//
// Pieces with a caliper width below D_min (A 393) are not created. Mass, momentum and angular
// momentum are conserved: pieces keep the orientation and the rigid-body velocity field of the parent.

namespace
{
    double polygon_centroid(const vector<double> &x, const vector<double> &y, double &cx, double &cy)
    {
        const int n=int(x.size());
        double A=0.0;
        cx=cy=0.0;
        for(int q=0;q<n;++q)
        {
        const int q2=(q+1)%n;
        const double cr = x[q]*y[q2]-x[q2]*y[q];
        A  += 0.5*cr;
        cx += (x[q]+x[q2])*cr;
        cy += (y[q]+y[q2])*cr;
        }
        if(fabs(A)>0.0)
        {
        cx/=(6.0*A);
        cy/=(6.0*A);
        }
        return A;
    }

    double caliper(const vector<double> &x, const vector<double> &y)
    {
        const int n=int(x.size());
        double w=1.0e20;
        for(int q=0;q<n;++q)
        {
        const int q2=(q+1)%n;
        const double ex=x[q2]-x[q], ey=y[q2]-y[q];
        const double len=sqrt(ex*ex+ey*ey);
        if(len<=0.0)
        continue;
        double dm=0.0;
        for(int r=0;r<n;++r)
        dm = MAX(dm, fabs(-ey*(x[r]-x[q]) + ex*(y[r]-y[q]))/len);
        w = MIN(w,dm);
        }
        return w>1.0e19 ? 0.0 : w;
    }
}

void fnpf_ice::clip_halfplane(const vector<double> &x, const vector<double> &y, double nx, double ny, double s, double sgn,
                              vector<double> &ox, vector<double> &oy)
{
    // keep sgn*(n.x - s) <= 0 (Sutherland-Hodgman, one plane)
    ox.clear();
    oy.clear();
    const int n=int(x.size());

    for(int q=0;q<n;++q)
    {
        const int q2=(q+1)%n;
        const double d0 = sgn*(nx*x[q]  + ny*y[q]  - s);
        const double d1 = sgn*(nx*x[q2] + ny*y[q2] - s);

        if(d0<=0.0)
        {
        ox.push_back(x[q]);
        oy.push_back(y[q]);
        }

        if((d0<0.0 && d1>0.0) || (d0>0.0 && d1<0.0))
        {
        const double t = d0/(d0-d1);
        ox.push_back(x[q] + t*(x[q2]-x[q]));
        oy.push_back(y[q] + t*(y[q2]-y[q]));
        }
    }
}

double fnpf_ice::chord(const vector<double> &x, const vector<double> &y, double nx, double ny, double s)
{
    // length of the intersection of the line n.x = s with a convex polygon
    const int n=int(x.size());
    double tmin=1.0e20, tmax=-1.0e20;
    const double tx=-ny, ty=nx;

    for(int q=0;q<n;++q)
    {
        const int q2=(q+1)%n;
        const double d0 = nx*x[q]  + ny*y[q]  - s;
        const double d1 = nx*x[q2] + ny*y[q2] - s;

        if((d0<=0.0 && d1>0.0) || (d0>0.0 && d1<=0.0))
        {
        const double t = d0/(d0-d1);
        const double xi = x[q] + t*(x[q2]-x[q]);
        const double yi = y[q] + t*(y[q2]-y[q]);
        const double u = tx*xi + ty*yi;
        tmin = MIN(tmin,u);
        tmax = MAX(tmax,u);
        }
    }
    return tmax>tmin ? tmax-tmin : 0.0;
}

void fnpf_ice::breaking_moments(lexer *p)
{
    // all ranks: cut moments from the local footprint cells, summed over the ranks
    const int ncut = ndir*noff;
    Mcut.assign(ncut*floe.size(),0.0);

    if(!(breakflag&1))
    return;

    // cut offsets along each direction, from the planform extent (identical on all ranks)
    vector<double> soff(ncut*floe.size(),0.0);

    for(size_t f=0; f<floe.size(); ++f)
    {
        const fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0)
        continue;

        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double nx = cos(th), ny = sin(th);
            double smin=1.0e20, smax=-1.0e20;

            for(size_t q=0; q<fl.wx.size(); ++q)
            {
            const double d = nx*(fl.wx[q]-fl.x[0]) + ny*(fl.wy[q]-fl.x[1]);
            smin = MIN(smin,d);
            smax = MAX(smax,d);
            }

            for(int o=0; o<noff; ++o)
            soff[ncut*f + k*noff + o] = smin + (smax-smin)*double(o+1)/double(noff+1);
        }
    }

    for(auto &e : cell)
    {
        const fnpf_ice_floe &fl = floe[e.f];

        if(fl.type!=0)
        continue;

        const double rx = e.xc - fl.x[0];
        const double ry = e.yc - fl.x[1];
        const double az = fl.acc[2] + fl.alp[0]*ry - fl.alp[1]*rx;
        const double qA = e.phi*e.area*(e.pl - fl.pw - fl.rho*fl.h*az);

        double *M = &Mcut[ncut*e.f];
        const double *so = &soff[ncut*e.f];

        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double d = cos(th)*rx + sin(th)*ry;

            for(int o=0; o<noff; ++o)
            {
            const double s = so[k*noff+o];
            if(d>s)
            M[k*noff+o] += qA*(d-s);
            }
        }
    }

    if(p->mpi_size>1 && !Mcut.empty())
    MPI_Allreduce(MPI_IN_PLACE, Mcut.data(), int(Mcut.size()), MPI_DOUBLE, MPI_SUM, comm);

    // keep the offsets for breaking_decide
    Mcut.insert(Mcut.end(), soff.begin(), soff.end());
}

void fnpf_ice::breaking_decide(lexer *p)
{
    // rank 0: pick at most one cut per floe
    splits.clear();

    const int ncut = ndir*noff;
    const size_t nf = floe.size();

    vector<split> best(nf);
    vector<double> ratio(nf,0.0);

    auto body_cut = [&](const fnpf_ice_floe &fl, double nx, double ny, double s, split &sp)
    {
        // world horizontal cut through the centre of mass frame -> body frame
        double bxn = fl.R[0][0]*nx + fl.R[1][0]*ny;
        double byn = fl.R[0][1]*nx + fl.R[1][1]*ny;
        const double len = sqrt(bxn*bxn + byn*byn);
        if(len<1.0e-12)
        return 0;
        sp.nbx = bxn/len;
        sp.nby = byn/len;
        sp.s = s/len;

        // both pieces at least D_min wide
        vector<double> ax,ay,bx2,by2;
        clip_halfplane(fl.bx,fl.by,sp.nbx,sp.nby,sp.s, 1.0,ax,ay);
        clip_halfplane(fl.bx,fl.by,sp.nbx,sp.nby,sp.s,-1.0,bx2,by2);
        if(ax.size()<3 || bx2.size()<3)
        return 0;
        if(is2D)
        {
        const double la = *max_element(ax.begin(),ax.end()) - *min_element(ax.begin(),ax.end());
        const double lb = *max_element(bx2.begin(),bx2.end()) - *min_element(bx2.begin(),bx2.end());
        return (la>=Dmin && lb>=Dmin) ? 1 : 0;
        }
        return (caliper(ax,ay)>=Dmin && caliper(bx2,by2)>=Dmin) ? 1 : 0;
    };

    // flexural
    if(breakflag&1)
    for(size_t f=0; f<nf; ++f)
    {
        const fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0 || fl.Awet<=0.0)
        continue;

        // planform at the centre of mass level, relative to it
        const int nv = int(fl.bx.size());
        vector<double> px(nv),py(nv);
        for(int q=0;q<nv;++q)
        {
        px[q] = fl.R[0][0]*fl.bx[q] + fl.R[0][1]*fl.by[q];
        py[q] = fl.R[1][0]*fl.bx[q] + fl.R[1][1]*fl.by[q];
        }

        for(int k=0; k<ndir; ++k)
        {
            const double th = PI*double(k)/double(ndir);
            const double nx = cos(th), ny = sin(th);

            for(int o=0; o<noff; ++o)
            {
                const double M = Mcut[ncut*f + k*noff + o];
                const double s = Mcut[ncut*nf + ncut*f + k*noff + o];
                const double L = is2D ? fl.width2D : chord(px,py,nx,ny,s);

                if(L<=0.0)
                continue;

                const double sig = 6.0*fabs(M)/(fl.h*fl.h*L);
                const double r = sig/MAX(sigf,1.0e-20);

                if(r>=1.0 && r>ratio[f])
                {
                split sp;
                if(body_cut(fl,nx,ny,s,sp))
                {
                sp.f = int(f);
                sp.mech = 1;
                sp.val = sig;
                sp.lim = sigf;
                best[f] = sp;
                ratio[f] = r;
                }
                }
            }
        }
    }

    // contact splitting, on the time-averaged pair normal force (contact_average)
    if((breakflag&2) && !is2D)
    for(const auto &kv : cforce)
    for(int side=0; side<2; ++side)
    {
        const int fid = side==0 ? kv.first.first : kv.first.second;
        if(fid<0 || fid>=int(nf))
        continue;
        
        const size_t f = size_t(fid);
        const fnpf_ice_floe &fl = floe[f];
        if(fl.type!=0)
        continue;
        
        const cavg &ca = kv.second;
        const double D = 2.0*sqrt(fl.area/PI);
        const double Fs = Csplit*KIC*fl.h*sqrt(D);
        const double r = ca.F/MAX(Fs,1.0e-20);
        
        if(r>=1.0 && r>ratio[f])
        {
        // cut along the contact normal through the contact point
        const double nx = -ca.ny, ny = ca.nx;
        const double s = nx*(ca.px-fl.x[0]) + ny*(ca.py-fl.x[1]);
        split sp;
        if(body_cut(fl,nx,ny,s,sp))
        {
        sp.f = int(f);
        sp.mech = 2;
        sp.val = ca.F;
        sp.lim = Fs;
        best[f] = sp;
        ratio[f] = r;
        }
        }
    }
    
    for(size_t f=0; f<nf; ++f)
    if(ratio[f]>=1.0)
    splits.push_back(best[f]);
}

void fnpf_ice::breaking_apply(lexer *p, ghostcell *pgc)
{
    int ns = (p->mpirank==0) ? int(splits.size()) : 0;

    if(p->mpi_size>1)
    {
    MPI_Bcast(&ns, 1, MPI_INT, 0, pgc->mpi_comm);

    vector<double> buf(7*ns);
    if(p->mpirank==0)
    for(int n=0;n<ns;++n)
    {
    const split &sp = splits[n];
    double *b=&buf[7*n];
    b[0]=sp.f; b[1]=sp.mech; b[2]=sp.nbx; b[3]=sp.nby; b[4]=sp.s; b[5]=sp.val; b[6]=sp.lim;
    }

    if(ns>0)
    MPI_Bcast(buf.data(), 7*ns, MPI_DOUBLE, 0, pgc->mpi_comm);

    if(p->mpirank>0)
    {
    splits.resize(ns);
    for(int n=0;n<ns;++n)
    {
    const double *b=&buf[7*n];
    splits[n].f=int(b[0]); splits[n].mech=int(b[1]); splits[n].nbx=b[2]; splits[n].nby=b[3];
    splits[n].s=b[4]; splits[n].val=b[5]; splits[n].lim=b[6];
    }
    }
    }

    const double t = p->simtime + p->dt;

    for(int n=0;n<ns;++n)
    {
        const split &sp = splits[n];
        const int nid = split_floe(p,size_t(sp.f),sp.nbx,sp.nby,sp.s,sp.mech);

        if(nid>=0)
        {
        ++nbreak;
        
        // the load that broke the floe is released
        const int pid = floe[sp.f].id;
        for(auto it=cforce.begin(); it!=cforce.end();)
        {
            if(it->first.first==pid || it->first.second==pid)
            it = cforce.erase(it);
            else
            ++it;
        }

        if(p->mpirank==0 && breakout.is_open())
        breakout<<setprecision(9)<<t<<" "<<floe[sp.f].id<<" "<<nid<<" "<<sp.mech<<" "<<sp.val<<" "<<sp.lim<<" "
                <<floe[sp.f].area<<" "<<floe.back().area<<endl;
        }
    }

    if(ns>0 && p->mpirank==0)
    cout<<"FNPF ice: "<<ns<<" floe(s) broken at t = "<<t<<", "<<nfloe<<" floes"<<endl;

    splits.clear();
}

int fnpf_ice::split_floe(lexer *p, size_t f, double nbx, double nby, double s, int mech)
{
    fnpf_ice_floe parent = floe[f];

    vector<double> ax,ay,bx,by;
    clip_halfplane(parent.bx,parent.by,nbx,nby,s, 1.0,ax,ay);
    clip_halfplane(parent.bx,parent.by,nbx,nby,s,-1.0,bx,by);

    if(ax.size()<3 || bx.size()<3)
    return -1;

    quat_to_matrix(parent.q,parent.R);

    auto make_piece = [&](const vector<double> &px, const vector<double> &py, fnpf_ice_floe &pc)
    {
        double cx,cy;
        polygon_centroid(px,py,cx,cy);

        pc = parent;
        pc.bx = px;
        pc.by = py;
        geometry(pc);   // re-centres on the piece centroid

        // world position and rigid-body velocity of the piece centre
        const double r[3] = {parent.R[0][0]*cx + parent.R[0][1]*cy,
                             parent.R[1][0]*cx + parent.R[1][1]*cy,
                             parent.R[2][0]*cx + parent.R[2][1]*cy};
        for(int a=0;a<3;++a)
        pc.x[a] = parent.x[a] + r[a];

        pc.v[0] = parent.v[0] + parent.w[1]*r[2] - parent.w[2]*r[1];
        pc.v[1] = parent.v[1] + parent.w[2]*r[0] - parent.w[0]*r[2];
        pc.v[2] = parent.v[2] + parent.w[0]*r[1] - parent.w[1]*r[0];

        pc.Fc[0]=pc.Fc[1]=0.0;
        pc.ncontact=0;
        kinematics(pc);
    };

    fnpf_ice_floe pa,pb;
    make_piece(ax,ay,pa);
    make_piece(bx,by,pb);

    pa.id = parent.id;
    pb.id = int(floe.size());

    floe[f] = pa;
    floe.push_back(pb);

    Yn.resize(floe.size());
    D1.resize(floe.size());
    D2.resize(floe.size());
    D3.resize(floe.size());
    ++nfloe;

    (void)mech;
    return pb.id;
}

void fnpf_ice::contact_average(lexer *p)
{
    // Rank 0, every step. The non-smooth contact gives impulses, an impact is resolved within one or two
    // steps and impulse/dt depends on dt and on where the impact falls within the step. For splitting
    // the pair normal force is averaged over t_c (A 396) with an exponential filter:
    //     F_avg <- F_avg + a*(F - F_avg),  a = min(1, dt/t_c)
    // An impact with impulse P gives F_avg ~ P/t_c, a sustained load gives F_avg = F.
    const double a = (tcon>0.0) ? MIN(1.0, p->dt/tcon) : 1.0;
    
    for(auto &kv : cforce)
    kv.second.F *= (1.0-a);
    
    for(const auto &rc : pcontact->records())
    {
        if(rc.a<0 || rc.a>=int(cmap.size()))
        continue;
        
        const int ida = floe[cmap[rc.a]].id;
        const int idb = (rc.b>=0 && rc.b<int(cmap.size())) ? floe[cmap[rc.b]].id : rc.b;
        
        cavg &ca = cforce[make_pair(ida,idb)];
        ca.F  += a*rc.Fn;
        ca.px = rc.px;
        ca.py = rc.py;
        ca.nx = rc.nx;
        ca.ny = rc.ny;
    }
    
    for(auto it=cforce.begin(); it!=cforce.end();)
    {
        if(it->second.F < 1.0e-6)
        it = cforce.erase(it);
        else
        ++it;
    }
}
