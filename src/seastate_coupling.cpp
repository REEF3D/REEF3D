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
Architect: Hans Bihs
--------------------------------------------------------------------*/


#include"seastate_coupling.h"
#include"seastate_remote.h"
#include"seastate_f.h"
#include"seastate_source.h"
#include"seastate_roller.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_dispersion.h"
#include"seastate_vtp.h"
#include"fdm_seastate.h"
#include"lexer.h"
#include"ghostcell.h"
#include"met_forcing.h"
#include<algorithm>
#include<iomanip>
#include<mpi.h>
#include<cmath>
#include<vector>

seastate_coupling::seastate_coupling(lexer *p) : pwave(nullptr), pnet(nullptr),
                Hsw(p),Tpw(p),Fw(p),Gw(p),divM(p),Sxx(p),Sxy(p),Syy(p),Mx(p),My(p),Fdx(p),Fdy(p),Q(p),B(p),
                tlast(0.0), tnext(0.0), tprint(0.0), nstep(0), remote(false), host_wind(false), t_wait(0.0)
{
}

bool seastate_coupling::surfbeat() const
{
    return pwave!=nullptr && pwave->surfbeat();
}

seastate_coupling::~seastate_coupling()
{
    delete pnet;
    delete pwave;
}

void seastate_coupling::ini_wave(lexer *p, ghostcell *pgc, const char *host)
{
    if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D::SEASTATE coupled with "<<host<<(remote ? " (A 750 2, its own grid): " : " (A 750 1): ");
    if(p->A770==1)
    cout<<"surfbeat (A 770 1), wave step every "<<host<<" step";
    else
    cout<<"coupling interval (wave time step) "<<p->A706<<" s";
    cout<<", wave forcing "<<(p->A751==1 ? "radiation stress" : "vortex force")<<", ramp-up "<<p->A752<<" s, feedback "
        <<(p->A753==0 ? "none" : (p->A753==1 ? "water level" : "water level and currents"))<<endl;
    }

    pwave = new seastate_f(p,pgc);

    if(remote)
    pwave->ini_remote(p,pgc);
    else
    pwave->ini_coupled(p,pgc);

    // wave force and mass flux in the SEASTATE VTP files (from the first wave step on)
    pwave->printer()->add_field("Fx",&Fw);
    pwave->printer()->add_field("Fy",&Gw);
    pwave->printer()->add_field("MStx",&Mx);
    pwave->printer()->add_field("MSty",&My);

    // net source terms for the vortex-force dissipation term (with the roller, the roller dissipation replaces the
    // breaking dissipation). The wind input enters when the host gets the full wind stress (met forcing W 200,
    // NHFLOW wind A 570): the momentum the wind gives to the waves is then taken off the mean flow (the
    // wave-supported stress) and comes back through the dissipation, else a growing wind sea counts it twice;
    // without a wind stress on the host it stays out (the waves hand over only what they lose)
    if(p->A751==2 && pwave->source()!=nullptr)
    {
    seastate_source_param prm = pwave->source()->param();
    prm.wind = prm.wind && (remote ? host_wind : (p->W200>0 || p->A570>0));

    if(pwave->roller()!=nullptr)
    prm.breaking = false;

    if(prm.any())
    pnet = new seastate_source(*pwave->e->grid,prm);
    }

    tlast = p->simtime;
    tnext = p->simtime + p->A706;
    tprint = p->simtime + (p->P182>0.0 ? p->P182 : 0.0);
    nstep = 0;
}

double seastate_coupling::ramp(lexer *p) const
{
    return (p->A752>0.0) ? std::min(1.0,std::max(p->simtime,0.0)/p->A752) : 1.0;
}

bool seastate_coupling::due(lexer *p, double &dt) const
{
    dt = p->simtime - tlast;

    if(surfbeat())
    return dt>0.0;

    return !(p->simtime < tnext - 1.0e-9*p->A706);
}

void seastate_coupling::wave_step(lexer *p, ghostcell *pgc, double dt)
{
    pwave->step_coupled(p,pgc,dt);

    tlast = p->simtime;
    tnext = tlast + p->A706;
    ++nstep;

    if(p->mpirank==0 && (!surfbeat() || p->count%p->P12==0))
    cout<<"SEASTATE wave step "<<nstep<<" at simtime "<<p->simtime<<" (dt "<<dt<<" s), forcing ramp "<<ramp(p)<<endl;
}

void seastate_coupling::print(lexer *p, ghostcell *pgc)
{
    if(p->P10<1)
    return;

    // coupling interval: every wave step
    if(!surfbeat())
    {
    pwave->printer()->print2D(p,pwave->e,pgc);
    return;
    }

    // surfbeat: own print interval P 182 [s] or P 181 [iterations]
    if(p->P182>0.0)
    {
        if(p->simtime>=tprint-1.0e-9)
        {
        pwave->printer()->print2D(p,pwave->e,pgc);
        tprint += p->P182;
        }
    }
    else if(p->P181>0 && p->count%p->P181==0)
    pwave->printer()->print2D(p,pwave->e,pgc);
}

void seastate_coupling::environment(lexer *p, slice &eta, slice &WL, slice &U, slice &V)
{
    fdm_seastate *e = pwave->e;

    IMALOOP
    JMALOOP
    {
    const bool inside = i+p->origin_i>=0 && i+p->origin_i<p->gknox && j+p->origin_j>=0 && j+p->origin_j<p->gknoy;

        if(!inside || e->wet0(i,j)==0)
        continue;

        // water level: total water depth of the host
        if(p->A753>=1)
        {
        e->eta(i,j)   = eta(i,j);
        e->depth(i,j) = std::max(WL(i,j),0.0);
        }

        // current: depth-averaged Eulerian velocity
        if(p->A753==2)
        {
        double u = U(i,j), v = V(i,j);

            // radiation stress: the host velocity is the mass transport velocity
            if(p->A751==1 && WL(i,j)>p->A705)
            {
            u -= Mx(i,j)/WL(i,j);
            v -= My(i,j)/WL(i,j);
            }

        e->U(i,j) = u;
        e->V(i,j) = v;
        }
    }
}

// derivatives between wave-active cells: central, one-sided next to dry cells and at the domain edge
double seastate_coupling::ddx(lexer *p, slice &f, int ii, int jj)
{
    fdm_seastate *e = pwave->e;
    i = ii;
    j = jj;

    const bool w = e->wet(i-1,j)==1, o = e->wet(i+1,j)==1;

    if(w && o)
    return (f(i+1,j)-f(i-1,j))/(p->XP[IP1]-p->XP[IM1]);

    if(o)
    return (f(i+1,j)-f(i,j))/(p->XP[IP1]-p->XP[IP]);

    if(w)
    return (f(i,j)-f(i-1,j))/(p->XP[IP]-p->XP[IM1]);

    return 0.0;
}

double seastate_coupling::ddy(lexer *p, slice &f, int ii, int jj)
{
    fdm_seastate *e = pwave->e;
    i = ii;
    j = jj;

    if(p->j_dir==0)
    return 0.0;

    const bool s = e->wet(i,j-1)==1, n = e->wet(i,j+1)==1;

    if(s && n)
    return (f(i,j+1)-f(i,j-1))/(p->YP[JP1]-p->YP[JM1]);

    if(n)
    return (f(i,j+1)-f(i,j))/(p->YP[JP1]-p->YP[JP]);

    if(s)
    return (f(i,j)-f(i,j-1))/(p->YP[JP]-p->YP[JM1]);

    return 0.0;
}

void seastate_coupling::forces(lexer *p, ghostcell *pgc, slice &U, slice &V)
{
    fdm_seastate *e = pwave->e;
    const seastate_grid &g = *e->grid;
    const seastate_roller *proll = pwave->roller();
    const double grav = seastate_gravity;
    const bool vf = (p->A751==2);

    std::vector<double> P(g.nbin), D(g.nbin);

    // spectral integrals per wave-active cell (incl. the ghost cells: the spectra are exchanged)
    IMALOOP
    JMALOOP
    {
    Sxx(i,j)=Sxy(i,j)=Syy(i,j)=Mx(i,j)=My(i,j)=Fdx(i,j)=Fdy(i,j)=Q(i,j)=B(i,j)=0.0;

        if(e->wet(i,j)==0 || e->N->spec(i,j)==nullptr)
        continue;

    const float *N  = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cc = e->cg->spec(i,j);
    const double d  = e->depth(i,j);

        if(vf && pnet!=nullptr)
        {
            // the wind of the cell (met forcing)
            if(pnet->param().wind && pwave->wind_speed()!=nullptr)
            pnet->set_wind((*pwave->wind_speed())(i,j),(*pwave->wind_direction())(i,j));

        pnet->compute(N,d,kc,cc,P.data(),D.data());
        }

        for(int l=0; l<g.nsig; ++l)
        {
        const double sig = g.sig[l];
        const double k = double(kc[l]);

            if(!(k>0.0))
            continue;

        const double n = double(cc[l])*k/sig;
        const double A = seastate_refraction(sig,k,d);          // sig/sinh(2kd)
        const double w = g.dsig[l]*g.dtheta;

            for(int m=0; m<g.ndir; ++m)
            {
            const int bb = g.bin(l,m);
            const double cs = g.costh[m], sn = g.sinth[m];
            const double E = sig*double(N[bb])*w*g.wth[m];      // [m^2]

            Mx(i,j) += grav*E*k/sig*cs;
            My(i,j) += grav*E*k/sig*sn;

                if(!vf)
                {
                Sxx(i,j) += grav*E*(n*(cs*cs+1.0)-0.5);
                Sxy(i,j) += grav*E*n*cs*sn;
                Syy(i,j) += grav*E*(n*(sn*sn+1.0)-0.5);
                }
                else
                {
                Q(i,j) += grav*E*(n-0.5);
                B(i,j) += grav*E/sig*k*A;

                    if(pnet!=nullptr)
                    {
                    const double S = (P[bb] - D[bb]*double(N[bb]))*w*g.wth[m];   // net source dN/dt dsig dtheta
                    Fdx(i,j) -= grav*S*k*cs;
                    Fdy(i,j) -= grav*S*k*sn;
                    }
                }
            }
        }

        // roller (surfbeat, one frequency): stress, mass flux 2 g R/c, dissipation force
        if(proll!=nullptr && proll->spec(i,j)!=nullptr)
        {
        const float *r = proll->spec(i,j);
        const double sig = g.sig[0];
        const double k = double(kc[0]);
        const double c = k>0.0 ? sig/k : 0.0;

            if(c>0.0)
            for(int m=0; m<g.ndir; ++m)
            {
            const double cs = g.costh[m], sn = g.sinth[m];
            const double Rm = double(r[m])*g.dtheta;             // [m^2]

            Mx(i,j) += 2.0*grav*Rm/c*cs;
            My(i,j) += 2.0*grav*Rm/c*sn;

                if(!vf)
                {
                Sxx(i,j) += grav*Rm*cs*cs;
                Sxy(i,j) += grav*Rm*cs*sn;
                Syy(i,j) += grav*Rm*sn*sn;
                }
                else
                {
                const double Dr = 2.0*grav*proll->beta/c*Rm;      // roller dissipation [m^2/s]
                Fdx(i,j) += grav*Dr/sig*k*cs;
                Fdy(i,j) += grav*Dr/sig*k*sn;
                }
            }
        }
    }

    // forces in the host cells
    SLICELOOP4
    {
    Fw(i,j)=Gw(i,j)=divM(i,j)=0.0;

        if(e->wet(i,j)==0)
        continue;

        if(!vf)
        {
        Fw(i,j) = -(ddx(p,Sxx,i,j) + ddy(p,Sxy,i,j));
        Gw(i,j) = -(ddx(p,Sxy,i,j) + ddy(p,Syy,i,j));
        }
        else
        {
        const int ic=i, jc=j;

        // div(M) in flux form: face values between wave-active cells, no Stokes transport
        // through the domain edges, walls or the shoreline, so the mass source integrates to zero
        const double fe = (e->wet(i+1,j)==1) ? 0.5*(Mx(i,j)+Mx(i+1,j)) : 0.0;
        const double fw = (e->wet(i-1,j)==1) ? 0.5*(Mx(i,j)+Mx(i-1,j)) : 0.0;
        const double fn = (p->j_dir==1 && e->wet(i,j+1)==1) ? 0.5*(My(i,j)+My(i,j+1)) : 0.0;
        const double fs = (p->j_dir==1 && e->wet(i,j-1)==1) ? 0.5*(My(i,j)+My(i,j-1)) : 0.0;
        const double dM  = (fe-fw)/p->DXN[IP] + (fn-fs)/p->DYN[JP];

        const double chi = ddx(p,V,ic,jc) - ddy(p,U,ic,jc);
        const double dQx = ddx(p,Q,ic,jc), dQy = ddy(p,Q,ic,jc);
        i = ic;
        j = jc;

        divM(i,j) = dM;
        Fw(i,j) = Fdx(i,j) + B(i,j)*e->ddx(i,j) - dQx + chi*My(i,j) - U(i,j)*dM;
        Gw(i,j) = Fdy(i,j) + B(i,j)*e->ddy(i,j) - dQy - chi*Mx(i,j) - V(i,j)*dM;
        }
    }

    // sea state on the host grid (the same grid in-process)
    IMALOOP
    JMALOOP
    {
    Hsw(i,j) = e->wet(i,j)==1 ? e->Hs(i,j) : 0.0;
    Tpw(i,j) = e->wet(i,j)==1 ? e->Tp(i,j) : 0.0;
    }

    if(pmetw!=nullptr)
    pmetw->waves_update(p);
}

bool seastate_coupling::longwave_bc(lexer *p, int jj, double h, double eta1, double v1, double &etag, double &ug, double &vg)
{
    if(pwave==nullptr || !pwave->surfbeat() || p->origin_i!=0 || !(h>p->A705))
    return false;

    double z=0.0, qx=0.0, qy=0.0;
    pwave->longwave(p,jj,p->simtime,z,qx,qy);

    etag = eta1;
    ug = qx/h - std::sqrt(std::fabs(p->W22)/h)*(eta1 - z);
    vg = v1;

    return true;
}

// ------------------------------------------------------------------ remote waves (A 750 2)

void seastate_coupling::centres(lexer *p, ghostcell *pgc, std::vector<double> &xs, std::vector<double> &ys)
{
    // global cell centres of this model (rectilinear grid): the largest value of every global index
    std::vector<double> lx(p->gknox,-1.0e300), ly(std::max(p->gknoy,1),-1.0e300);

    for(i=0; i<p->knox; ++i)
    lx[i+p->origin_i] = p->XP[IP];

    for(j=0; j<p->knoy; ++j)
    ly[j+p->origin_j] = p->YP[JP];

    xs.assign(lx.size(),0.0);
    ys.assign(ly.size(),0.0);
    MPI_Allreduce(lx.data(),xs.data(),int(lx.size()),MPI_DOUBLE,MPI_MAX,pgc->mpi_comm);
    MPI_Allreduce(ly.data(),ys.data(),int(ly.size()),MPI_DOUBLE,MPI_MAX,pgc->mpi_comm);
}

void seastate_coupling::gather(lexer *p, ghostcell *pgc, const std::vector<slice*> &f, std::vector<double> &out, double flag)
{
    // the fields of the own cells into global arrays [field][I*ny + J] behind the header (flag, simtime), on rank 0
    const int ny = std::max(p->gknoy,1);
    const size_t N = size_t(p->gknox)*size_t(ny);
    std::vector<double> loc(2+f.size()*N,0.0);

    for(size_t q=0; q<f.size(); ++q)
    SLICELOOP4
    loc[2+q*N+size_t(i+p->origin_i)*ny+size_t(j+p->origin_j)] = (*f[q])(i,j);

    out.assign(loc.size(),0.0);
    MPI_Reduce(loc.data(),out.data(),int(loc.size()),MPI_DOUBLE,MPI_SUM,0,pgc->mpi_comm);

    out[0] = flag;
    out[1] = p->simtime;
}

void seastate_coupling::exchange(lexer *p, ghostcell *pgc, std::vector<double> &out, std::vector<double> &in, int tag_out, int tag_in)
{
    // root to root over the intercommunicator (send and receive at once: no order between the apps), then to all ranks
    if(p->mpirank==0)
    {
    const double t0 = MPI_Wtime();
    MPI_Request rq;
    MPI_Isend(out.data(),int(out.size()),MPI_DOUBLE,0,tag_out,seastate_remote::ic,&rq);
    MPI_Recv(in.data(),int(in.size()),MPI_DOUBLE,0,tag_in,seastate_remote::ic,MPI_STATUS_IGNORE);
    MPI_Wait(&rq,MPI_STATUS_IGNORE);
    t_wait += MPI_Wtime() - t0;
    }

    MPI_Bcast(in.data(),int(in.size()),MPI_DOUBLE,0,pgc->mpi_comm);
}

void seastate_coupling::bilinear(const std::vector<double> &xs, const std::vector<double> &ys, double x, double y,
                                 int &a, int &b, double &wx, double &wy) const
{
    // cell (a, b) of the lattice of centres and the weights of a+1, b+1 (clamped: constant beyond the outer centres)
    const int nx = int(xs.size()), ny = int(ys.size());

    a = 0;
    wx = 0.0;

    if(nx>1)
    {
    a = int(std::upper_bound(xs.begin(),xs.end(),x) - xs.begin()) - 1;
    a = std::min(std::max(a,0),nx-2);
    wx = std::min(std::max((x - xs[a])/(xs[a+1] - xs[a]),0.0),1.0);
    }

    b = 0;
    wy = 0.0;

    if(ny>1)
    {
    b = int(std::upper_bound(ys.begin(),ys.end(),y) - ys.begin()) - 1;
    b = std::min(std::max(b,0),ny-2);
    wy = std::min(std::max((y - ys[b])/(ys[b+1] - ys[b]),0.0),1.0);
    }
}

void seastate_coupling::remote_ini(lexer *p, ghostcell *pgc, slice &eta, slice &WL, slice &U, slice &V)
{
    remote = true;

    // host-side checks (the waves check their keys in seastate_f)
    const char *msg = nullptr;

    if(p->A770==1)
    msg = "A 750 2: the surfbeat model (A 770 1) runs on the host grid and step (A 750 1)";
    else if(p->A751!=1 && p->A751!=2)
    msg = "A 751: the wave forcing must be 1 (radiation stress) or 2 (vortex force)";
    else if(p->A753<0 || p->A753>2)
    msg = "A 753: the feedback must be 0 (none), 1 (water level) or 2 (water level and currents)";
    else if(!(p->A706>0.0))
    msg = "A 706: the coupling interval must be positive";

    if(msg!=nullptr)
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<msg<<endl<<endl;

    pgc->final(true);
    }

    if(seastate_remote::ic==MPI_COMM_NULL)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! A 750 2: the remote waves need a SEASTATE model (A 10 7, A 750 2) in the same MPI job, e.g."<<endl
            <<"    mpirun -np 4 -wdir surge REEF3D : -np 2 -wdir waves REEF3D !!!"<<endl<<endl;

    pgc->final(true);
    }

    // set-up: the keys both sides must share, the grid of this model; then the grid of the waves
    std::vector<double> xs, ys;
    centres(p,pgc,xs,ys);

    const double mine[9] = {double(xs.size()), double(ys.size()), p->A706, double(p->A751), double(p->A753),
                            (p->W200>0 || p->A570>0) ? 1.0 : 0.0, p->W203, p->A780, double(p->A10)};
    double theirs[9] = {0,0,0,0,0,0,0,0,0};

    if(p->mpirank==0)
    {
    MPI_Sendrecv(mine,9,MPI_DOUBLE,0,4951,theirs,9,MPI_DOUBLE,0,4952,seastate_remote::ic,MPI_STATUS_IGNORE);

    std::vector<double> g(xs);
    g.insert(g.end(),ys.begin(),ys.end());
    rx.assign(size_t(theirs[0]),0.0);
    ry.assign(size_t(theirs[1]),0.0);
    std::vector<double> h(rx.size()+ry.size());
    MPI_Sendrecv(g.data(),int(g.size()),MPI_DOUBLE,0,4953,h.data(),int(h.size()),MPI_DOUBLE,0,4954,seastate_remote::ic,MPI_STATUS_IGNORE);
    std::copy(h.begin(),h.begin()+rx.size(),rx.begin());
    std::copy(h.begin()+rx.size(),h.end(),ry.begin());
    }

    MPI_Bcast(theirs,9,MPI_DOUBLE,0,pgc->mpi_comm);
    rx.resize(size_t(theirs[0]));
    ry.resize(size_t(theirs[1]));
    MPI_Bcast(rx.data(),int(rx.size()),MPI_DOUBLE,0,pgc->mpi_comm);
    MPI_Bcast(ry.data(),int(ry.size()),MPI_DOUBLE,0,pgc->mpi_comm);

    if(theirs[2]!=p->A706 || int(theirs[3])!=p->A751 || int(theirs[4])!=p->A753)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! A 750 2: host and waves need the same A 706, A 751 and A 753 (host "<<p->A706<<" "<<p->A751<<" "<<p->A753
            <<", waves "<<theirs[2]<<" "<<theirs[3]<<" "<<theirs[4]<<") !!!"<<endl<<endl;

    MPI_Abort(MPI_COMM_WORLD,1);
    }

    if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D::SEASTATE remote (A 750 2): waves on their own grid "<<rx.size()<<" x "<<ry.size()<<" cells ("
        <<(rx.size()>1 ? rx[1]-rx[0] : 0.0)<<" m), host "<<xs.size()<<" x "<<ys.size()<<"; exchange every "<<p->A706
        <<" s, wave forcing "<<(p->A751==1 ? "radiation stress" : "vortex force")<<", ramp-up "<<p->A752<<" s"<<endl;

        // the clock of each model: A 780, else W 203 (0: none)
        const double th = p->A780>0.0 ? p->A780 : p->W203, tw = theirs[7]>0.0 ? theirs[7] : theirs[6];

        if(th>0.0 && tw>0.0 && th!=tw)
        {
        std::ios fmt(nullptr);
        fmt.copyfmt(cout);
        cout<<"!!! A 750 2: the waves start at "<<std::fixed<<std::setprecision(6)<<tw<<", the host at "<<th
            <<" (A 780, else W 203): both models must start at the same date-time !!!"<<endl;
        cout.copyfmt(fmt);
        }
    }

    // buffers: host fields eta, WL, U, V out; wave force, div(M), Hs, Tp, wet in
    const size_t Nw = rx.size()*ry.size();
    rbuf.assign(2+6*Nw,0.0);

    tlast = p->simtime;
    tnext = p->simtime;
    nstep = 0;

    remote_step(p,pgc,eta,WL,U,V,true);
}

void seastate_coupling::remote_step(lexer *p, ghostcell *pgc, slice &eta, slice &WL, slice &U, slice &V, bool force)
{
    if(!force && p->simtime < tnext - 1.0e-9*p->A706)
    return;

    gather(p,pgc,{&eta,&WL,&U,&V},sbuf,1.0);
    exchange(p,pgc,sbuf,rbuf,4961,4962);

    // the wave fields at the host cells: bilinear between the wet wave cell centres
    const int nxw = int(rx.size()), nyw = int(ry.size());
    const size_t Nw = size_t(nxw)*size_t(nyw);
    const double *F[5] = {&rbuf[2], &rbuf[2+Nw], &rbuf[2+2*Nw], &rbuf[2+3*Nw], &rbuf[2+4*Nw]};
    const double *W = &rbuf[2+5*Nw];
    slice4 *T[5] = {&Fw, &Gw, &divM, &Hsw, &Tpw};

    IMALOOP
    JMALOOP
    {
    int a, b;
    double wx, wy, v[5] = {0,0,0,0,0}, ws = 0.0;

        // host cells outside the wave grid: no waves
        if(!covers(rx,ry,p->XP[IP],p->YP[JP]))
        {
            for(int n=0; n<5; ++n)
            (*T[n])(i,j) = 0.0;
        continue;
        }

    bilinear(rx,ry,p->XP[IP],p->YP[JP],a,b,wx,wy);

        for(int da=0; da<2; ++da)
        for(int db=0; db<2; ++db)
        {
        const double w = (da==0 ? 1.0-wx : wx)*(db==0 ? 1.0-wy : wy);
        const int ia = std::min(a+da,nxw-1), jb = std::min(b+db,nyw-1);
        const size_t q = size_t(ia)*nyw + jb;

            if(w<=0.0 || W[q]<0.5)
            continue;

            for(int n=0; n<5; ++n)
            v[n] += w*F[n][q];

        ws += w;
        }

        for(int n=0; n<5; ++n)
        (*T[n])(i,j) = ws>0.0 ? v[n]/ws : 0.0;
    }

    if(pmetw!=nullptr)
    pmetw->waves_update(p);

    tlast = p->simtime;
    tnext = p->simtime + p->A706;
    ++nstep;

    if(p->mpirank==0)
    {
    std::ios fmt(nullptr);
    fmt.copyfmt(cout);
    cout<<"SEASTATE remote exchange "<<nstep<<" at simtime "<<std::fixed<<std::setprecision(1)<<p->simtime<<" s (wave time "<<rbuf[1]
        <<" s), forcing ramp "<<std::setprecision(3)<<ramp(p)<<endl;
    cout.copyfmt(fmt);
    }
}

void seastate_coupling::remote_stop(lexer *p)
{
    // rank 0 only (from ghostcell::final): the stop flag in an environment message, the pending wave fields drained
    if(!remote || p->mpirank!=0 || seastate_remote::ic==MPI_COMM_NULL || sbuf.empty())
    return;

    sbuf[0] = 0.0;
    MPI_Request rq;
    MPI_Isend(sbuf.data(),int(sbuf.size()),MPI_DOUBLE,0,4961,seastate_remote::ic,&rq);
    MPI_Recv(rbuf.data(),int(rbuf.size()),MPI_DOUBLE,0,4962,seastate_remote::ic,MPI_STATUS_IGNORE);
    MPI_Wait(&rq,MPI_STATUS_IGNORE);
    remote = false;

    cout<<"SEASTATE remote: "<<nstep<<" exchanges, the host waited "<<t_wait<<" s for the waves"<<endl;
}

bool seastate_coupling::covers(const std::vector<double> &xs, const std::vector<double> &ys, double x, double y) const
{
    // the point lies inside the grid of these cell centres (half a cell beyond the outer centres)
    const int nx = int(xs.size()), ny = int(ys.size());
    const double hx0 = nx>1 ? 0.5*(xs[1]-xs[0]) : 0.0, hx1 = nx>1 ? 0.5*(xs[nx-1]-xs[nx-2]) : 0.0;
    const double hy0 = ny>1 ? 0.5*(ys[1]-ys[0]) : 0.0, hy1 = ny>1 ? 0.5*(ys[ny-1]-ys[ny-2]) : 0.0;

    if(nx>1 && (x < xs[0]-hx0 || x > xs[nx-1]+hx1))
    return false;

    if(ny>1 && (y < ys[0]-hy0 || y > ys[ny-1]+hy1))
    return false;

    return true;
}
