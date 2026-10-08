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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#include"acoustics_fwh.h"
#include"acoustics_fwh_kernel.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>
#include<climits>
#include<cstdio>
#include<iomanip>
#include<iostream>
#include<sys/stat.h>
#include<sys/types.h>

acoustics_fwh::acoustics_fwh(lexer *p, fdm*, ghostcell*) : pfwh(nullptr), pout(nullptr), zfs(0.0), rho0(p->W1), tprev(-1.0e20), nobs(p->U30)
{
    box[0] = p->U20_xs;
    box[1] = p->U20_xe;
    box[2] = p->U20_ys;
    box[3] = p->U20_ye;
    box[4] = p->U20_zs;
    box[5] = p->U20_ze;
}

acoustics_fwh::~acoustics_fwh()
{
    delete pfwh;
    delete[] pout;
}

void acoustics_fwh::start(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->simtime<p->U32)
    return;
    
    if(pfwh==nullptr)
    {
        // default observer time step: the first source time step (p->dt is already the next one)
        if(p->U31<=0.0 && tprev<-1.0e19)
        {
            tprev = p->simtime;
            return;
        }
        
        ini(p,a,pgc);
    }
    
    sample(p,a);
    pfwh->step(p->simtime,pv.data(),uv.data());
    write(p,pgc);
    box_integrals(p,a,pgc);
}

void acoustics_fwh::ini(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->j_dir==0)
    {
        if(p->mpirank==0)
        std::cout<<"U 10: the FW-H acoustic analogy needs a 3D grid"<<std::endl;
        exit(1);
    }
    
    // still water level is known only after the initialisation of phi
    zfs = p->U41>-1.0e19 ? p->U41 : p->phimean;
    const double dto = p->U31>0.0 ? p->U31 : p->simtime-tprev;
    
    // snap the box to the nearest grid nodes
    const double *N[3] = {p->XN,p->YN,p->ZN};
    const int K[3] = {p->knox,p->knoy,p->knoz};
    
    for(int q=0; q<6; ++q)
    box[q] = snap(pgc,N[q/2],K[q/2],box[q]);
    
    const double gmin[3] = {p->global_xmin,p->global_ymin,p->global_zmin};
    const double gmax[3] = {p->global_xmax,p->global_ymax,p->global_zmax};
    
    for(int q=0; q<3; ++q)
    if(!(box[2*q]<box[2*q+1]) || box[2*q]<=gmin[q] || box[2*q+1]>=gmax[q])
    {
        if(p->mpirank==0)
        std::cout<<"U 20: the FW-H box must have a positive size and lie inside the domain"<<std::endl;
        exit(1);
    }
    
    // U 22 end caps, snapped to grid nodes, from the face U 21 inwards
    capx.clear();
    
    if(p->U21>0 && p->U22>0)
    {
        const int od = (p->U21-1)/2;
        const double sgn = (p->U21-1)%2==1 ? -1.0 : 1.0;
        
        for(int k=0; k<p->U22; ++k)
        {
            const double xk = box[p->U21-1] + sgn*(p->U22>1 ? p->U22_d*double(k)/double(p->U22-1) : 0.0);
            capx.push_back(snap(pgc,N[od],K[od],xk));
        }
        
        if(capx.back()<=box[2*od] || capx.back()>=box[2*od+1])
        {
            if(p->mpirank==0)
            std::cout<<"U 22: the end caps must lie inside the box U 20"<<std::endl;
            exit(1);
        }
    }
    
    box_faces(p,a,pgc);
    
    pfwh = new fwh_permeable(p->U11,rho0,dto);
    
    for(const face &f : fc)
    {
        fwh_panel P;
        for(int q=0; q<3; ++q)
        {
            P.x[q] = f.x[q];
            P.n[q] = 0.0;
        }
        P.n[f.dir] = f.sgn;
        P.dS = f.dS;
        pfwh->add_panel(P);
    }
    
    pv.assign(fc.size(),0.0);
    uv.assign(3*fc.size(),0.0);
    
    for(int n=0; n<nobs; ++n)
    {
        const double xo[3] = {p->U30_x[n],p->U30_y[n],p->U30_z[n]};
        const int o = pfwh->add_observer(xo);
        
        if(p->U40>0)
        {
            const double xi[3] = {xo[0],xo[1],2.0*zfs-xo[2]};
            pfwh->add_image(o,xi,-1.0);
        }
    }
    
    dmin.resize(size_t(nobs));
    dmax.resize(size_t(nobs));
    kwrite.assign(size_t(nobs),LONG_MIN);
    
    for(int n=0; n<nobs; ++n)
    {
        double d0, d1;
        pfwh->delays(n,d0,d1);
        dmin[size_t(n)] = pgc->globalmin(d0);
        dmax[size_t(n)] = pgc->globalmax(d1);
    }
    
    const int npan = pgc->globalisum(int(fc.size()));
    
    if(p->mpirank==0)
    {
        std::cout<<"FW-H: "<<npan<<" panels, box "<<box[0]<<" "<<box[1]<<" "<<box[2]<<" "<<box[3]<<" "<<box[4]<<" "<<box[5]
                 <<", open face "<<p->U21<<", end caps "<<capx.size()<<", observer dt "<<dto<<", mirror "<<p->U40<<" at z = "<<zfs<<std::endl;
        
        mkdir("./REEF3D_CFD_Acoustics",0777);
        
        pout = new std::ofstream[size_t(nobs)];
        
        for(int n=0; n<nobs; ++n)
        {
            char name[200];
            snprintf(name,sizeof(name),"./REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-%i.dat",n+1);
            std::ofstream &out = pout[n];
            out.open(name);
            
            out<<"FW-H permeable surface, observer "<<n+1<<std::endl<<std::endl;
            out<<"x_coord     y_coord     z_coord"<<std::endl;
            out<<p->U30_x[n]<<" \t "<<p->U30_y[n]<<" \t "<<p->U30_z[n]<<std::endl<<std::endl;
            out<<"box: "<<box[0]<<" "<<box[1]<<" "<<box[2]<<" "<<box[3]<<" "<<box[4]<<" "<<box[5]<<", panels: "<<npan<<", open face: "<<p->U21<<std::endl;
            out<<"c0: "<<p->U11<<", rho0: "<<rho0<<", mirror: "<<p->U40<<" at z = "<<zfs<<std::endl;
            out<<"delay min/max: "<<dmin[size_t(n)]<<" "<<dmax[size_t(n)]<<std::endl<<std::endl;
            out<<"time \t p'"<<std::endl;
        }
        
        bout.open("./REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat");
        bout<<"FW-H box integrals: Lp_i = int p' n_i dS, Lm_i = int rho u_i u_n dS, V_i = int tau_ij n_j dS, M_i = int x_i rho u_n dS"<<std::endl;
        bout<<"force of the fluid on bodies at rest inside the box: F = -(Lp + Lm - V + dM/dt)"<<std::endl<<std::endl;
        bout<<"time \t Lpx \t Lpy \t Lpz \t Lmx \t Lmy \t Lmz \t Vx \t Vy \t Vz \t Mx \t My \t Mz"<<std::endl;
    }
}

void acoustics_fwh::box_integrals(lexer *p, fdm *a, ghostcell *pgc)
{
    const double c2 = p->U11*p->U11;
    const double *N[3] = {p->XN,p->YN,p->ZN};
    const double *P[3] = {p->XP,p->YP,p->ZP};
    double s[12] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    
    // staggered velocity component q at index c
    auto vel = [&](int q, const int *c) {return q==0 ? a->u(c[0],c[1],c[2]) : q==1 ? a->v(c[0],c[1],c[2]) : a->w(c[0],c[1],c[2]);};
    
    for(size_t n=0; n<fc.size(); ++n)
    {
        const face &f = fc[n];
        const int d = f.dir;
        const double rho = rho0 + pv[n]/c2;
        const double un  = f.sgn*uv[3*n+size_t(d)];
        
        const int lo[3] = {f.ii,f.jj,f.kk};
        int hi[3] = {f.ii,f.jj,f.kk};
        ++hi[d];
        
        const double mu = rho0*0.5*(a->visc(lo[0],lo[1],lo[2]) + a->visc(hi[0],hi[1],hi[2])
                                  + a->eddyv(lo[0],lo[1],lo[2]) + a->eddyv(hi[0],hi[1],hi[2]));
        
        for(int q=0; q<3; ++q)
        {
            s[q]   += (q==d ? f.sgn*pv[n] : 0.0)*f.dS;
            s[3+q] += rho*uv[3*n+size_t(q)]*un*f.dS;
            s[9+q] += f.x[q]*rho*un*f.dS;
            
            // viscous traction tau_qd n_d on the face
            double tau;
            
            if(q==d)
            {
                int m[3] = {lo[0],lo[1],lo[2]};
                --m[d];
                const double dudx = (vel(d,hi) - vel(d,m))/(N[d][lo[d]+2+marge] - N[d][lo[d]+marge]);
                tau = 2.0*mu*dudx;
            }
            else
            {
                // d u_q / d x_d: u_q at the two cell centres across the face
                int lom[3] = {lo[0],lo[1],lo[2]}, him[3] = {hi[0],hi[1],hi[2]};
                --lom[q];
                --him[q];
                const double uql = 0.5*(vel(q,lo) + vel(q,lom));
                const double uqh = 0.5*(vel(q,hi) + vel(q,him));
                const double duq = (uqh - uql)/(P[d][lo[d]+1+marge] - P[d][lo[d]+marge]);
                
                // d u_d / d x_q: u_d on the face, central along q
                int qp[3] = {lo[0],lo[1],lo[2]}, qm[3] = {lo[0],lo[1],lo[2]};
                ++qp[q];
                --qm[q];
                const double dud = (vel(d,qp) - vel(d,qm))/(P[q][lo[q]+1+marge] - P[q][lo[q]-1+marge]);
                
                tau = mu*(duq + dud);
            }
            
            s[6+q] += tau*f.sgn*f.dS;
        }
    }
    
    pgc->globalsum(s,12);
    
    if(p->mpirank==0)
    {
        bout<<std::setprecision(10)<<p->simtime;
        for(int q=0; q<12; ++q)
        bout<<" \t "<<s[q];
        bout<<std::endl;
    }
}

double acoustics_fwh::snap(ghostcell *pgc, const double *N, int K, double x)
{
    // nearest node of this rank, then the nearest over all ranks
    double best = N[marge], dist = fabs(N[marge]-x);
    
    for(int q=1; q<=K; ++q)
    if(fabs(N[q+marge]-x)<dist)
    {
        best = N[q+marge];
        dist = fabs(N[q+marge]-x);
    }
    
    const double gdist = pgc->globalmin(dist);
    
    return pgc->globalmax(dist==gdist ? best : -1.0e300);
}

void acoustics_fwh::box_faces(lexer *p, fdm *a, ghostcell *pgc)
{
    const double *N[3] = {p->XN,p->YN,p->ZN};
    const double *P[3] = {p->XP,p->YP,p->ZP};
    const double *D[3] = {p->DXN,p->DYN,p->DZN};
    const int K[3] = {p->knox,p->knoy,p->knoz};
    
    const double tol = 1.0e-9*(fabs(p->global_xmax-p->global_xmin) + fabs(p->global_ymax-p->global_ymin) + fabs(p->global_zmax-p->global_zmin));
    
    auto flagidx = [&](int ci, int cj, int ck) {return (ci-p->imin)*p->jmax*p->kmax + (cj-p->jmin)*p->kmax + ck-p->kmin;};
    
    int nbad=0;
    
    fc.clear();
    
    // planes of the surface: the box faces, or for the face U 21 the end caps
    struct plane {int dir; double X, sgn, w;};
    std::vector<plane> pl;
    const int od = p->U21>0 ? (p->U21-1)/2 : -1;
    const int os = p->U21>0 ? (p->U21-1)%2 : -1;
    const double nc = double(capx.size());
    
    for(int dir=0; dir<3; ++dir)
    for(int side=0; side<2; ++side)
    {
        const double sgn = side==0 ? -1.0 : 1.0;
        
        if(p->U21==2*dir+side+1)
        {
            for(double xk : capx)
            pl.push_back({dir,xk,sgn,1.0/nc});
        }
        else
        pl.push_back({dir,box[2*dir+side],sgn,1.0});
    }
    
    for(const plane &pn : pl)
    {
        const int dir = pn.dir;
        const double X = pn.X;
        
        // the face belongs to the rank with the node in its range [0,K)
        int nd=-1;
        for(int q=0; q<K[dir]; ++q)
        if(fabs(N[dir][q+marge]-X)<tol)
        nd=q;
        
        if(nd<0)
        continue;
        
        const int b = (dir+1)%3, c = (dir+2)%3;
        
        for(int tb=0; tb<K[b]; ++tb)
        for(int tc=0; tc<K[c]; ++tc)
        {
            if(N[b][tb+marge]<box[2*b]-tol || N[b][tb+1+marge]>box[2*b+1]+tol)
            continue;
            
            if(N[c][tc+marge]<box[2*c]-tol || N[c][tc+1+marge]>box[2*c+1]+tol)
            continue;
            
            // side faces with end caps: fraction of the closed surfaces that contain the panel
            double w = pn.w;
            
            if(nc>0.0 && dir!=od)
            {
                const int t = od==b ? tb : tc;
                const double lo = N[od][t+marge], hi = N[od][t+1+marge];
                int cnt=0;
                
                for(double xk : capx)
                if(os==1 ? hi<=xk+tol : lo>=xk-tol)
                ++cnt;
                
                w = double(cnt)/nc;
            }
            
            if(w<=0.0)
            continue;
            
            int idx[3];
            idx[dir] = nd-1;
            idx[b] = tb;
            idx[c] = tc;
            
            face f;
            f.ii = idx[0];
            f.jj = idx[1];
            f.kk = idx[2];
            f.dir = dir;
            f.sgn = pn.sgn;
            f.x[dir] = X;
            f.x[b] = P[b][tb+marge];
            f.x[c] = P[c][tc+marge];
            f.dS = w*D[b][tb+marge]*D[c][tc+marge];
            fc.push_back(f);
            
            // both cells must be fluid: no solid, no floating body, no air
            int up[3] = {idx[0],idx[1],idx[2]};
            ++up[dir];
            
            const bool solid = p->flag4[flagidx(idx[0],idx[1],idx[2])]<0 || p->flag4[flagidx(up[0],up[1],up[2])]<0;
            const bool body  = p->X10>0 && (a->fb(idx[0],idx[1],idx[2])<0.0 || a->fb(up[0],up[1],up[2])<0.0);
            const bool air   = a->phi(idx[0],idx[1],idx[2])<0.0 || a->phi(up[0],up[1],up[2])<0.0;
            
            if(solid || body || air)
            ++nbad;
        }
    }
    
    nbad = pgc->globalisum(nbad);
    
    if(p->mpirank==0 && nbad>0)
    std::cout<<"FW-H WARNING: "<<nbad<<" panels of the box U 20 touch a solid, a floating body or air; the box should lie in the water and enclose the bodies"<<std::endl;
}

void acoustics_fwh::sample(lexer *p, fdm *a)
{
    for(size_t n=0; n<fc.size(); ++n)
    {
        const face &f = fc[n];
        const int ci=f.ii, cj=f.jj, ck=f.kk;
        double pr, ux, uy, uz;
        
        if(f.dir==0)
        {
            pr = 0.5*(a->press(ci,cj,ck) + a->press(ci+1,cj,ck));
            ux = a->u(ci,cj,ck);
            uy = 0.25*(a->v(ci,cj,ck) + a->v(ci,cj-1,ck) + a->v(ci+1,cj,ck) + a->v(ci+1,cj-1,ck));
            uz = 0.25*(a->w(ci,cj,ck) + a->w(ci,cj,ck-1) + a->w(ci+1,cj,ck) + a->w(ci+1,cj,ck-1));
        }
        else if(f.dir==1)
        {
            pr = 0.5*(a->press(ci,cj,ck) + a->press(ci,cj+1,ck));
            ux = 0.25*(a->u(ci,cj,ck) + a->u(ci-1,cj,ck) + a->u(ci,cj+1,ck) + a->u(ci-1,cj+1,ck));
            uy = a->v(ci,cj,ck);
            uz = 0.25*(a->w(ci,cj,ck) + a->w(ci,cj,ck-1) + a->w(ci,cj+1,ck) + a->w(ci,cj+1,ck-1));
        }
        else
        {
            pr = 0.5*(a->press(ci,cj,ck) + a->press(ci,cj,ck+1));
            ux = 0.25*(a->u(ci,cj,ck) + a->u(ci-1,cj,ck) + a->u(ci,cj,ck+1) + a->u(ci-1,cj,ck+1));
            uy = 0.25*(a->v(ci,cj,ck) + a->v(ci,cj-1,ck) + a->v(ci,cj,ck+1) + a->v(ci,cj-1,ck+1));
            uz = a->w(ci,cj,ck);
        }
        
        // without the hydrostatic part rho0 g.(x - x_fs)
        const double phs = rho0*(p->W20*f.x[0] + p->W21*f.x[1] + p->W22*(f.x[2]-zfs));
        
        pv[n] = pr - p->pressgage - phs;
        uv[3*n]   = ux;
        uv[3*n+1] = uy;
        uv[3*n+2] = uz;
    }
}

void acoustics_fwh::write(lexer *p, ghostcell *pgc)
{
    // observer samples that became complete since the last call; the ranges come from global
    // quantities (source times, delay bounds over all ranks), so all ranks pack the same buffer
    const double dto = pfwh->dt_obs();
    const long k0 = pfwh->first_index();
    
    std::vector<long> ka(static_cast<size_t>(nobs),0), kb(static_cast<size_t>(nobs),-1);
    size_t ntot=0;
    
    for(int n=0; n<nobs; ++n)
    {
        const size_t o = size_t(n);
        double t0, t1;
        ka[o]=0;
        kb[o]=-1;
        
        if(!pfwh->complete_range(dmin[o],dmax[o],t0,t1))
        continue;
        
        kwrite[o] = std::max(kwrite[o],long(ceil(t0/dto)));
        ka[o] = kwrite[o];
        kb[o] = long(floor(t1/dto));
        
        if(kb[o]>=ka[o])
        ntot += size_t(kb[o]-ka[o]+1);
    }
    
    if(ntot==0)
    return;
    
    buf.assign(ntot,0.0);
    
    size_t m=0;
    for(int n=0; n<nobs; ++n)
    {
        const size_t o = size_t(n);
        const std::vector<double> &sig = pfwh->data(n);
        
        for(long k=ka[o]; k<=kb[o]; ++k)
        {
            const long ks = k-k0;
            buf[m] = ks>=0 && ks<long(sig.size()) ? sig[size_t(ks)] : 0.0;
            ++m;
        }
    }
    
    pgc->globalsum(buf.data(),int(ntot));
    
    m=0;
    for(int n=0; n<nobs; ++n)
    {
        const size_t o = size_t(n);
        
        for(long k=ka[o]; k<=kb[o]; ++k)
        {
            if(p->mpirank==0)
            pout[n]<<std::setprecision(10)<<double(k)*dto<<" \t "<<buf[m]<<"\n";
            ++m;
        }
        
        if(kb[o]>=ka[o])
        kwrite[o] = kb[o]+1;
        
        if(p->mpirank==0)
        pout[n].flush();
    }
}
