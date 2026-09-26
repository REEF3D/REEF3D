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

#include"dem_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4a.h"
#include<mpi.h>

// ---------------------------------------------------------------------------------------------
// fluid data at the particles: velocity, density, viscosity, voidage at the centroid,
// buoyancy from the volume quadrature (unresolved particles)
// ---------------------------------------------------------------------------------------------

void dem_f::fluid_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    const int nv = 14;
    vector<double> buf(nb*nv,0.0);
    dem_vec g(p->W20,p->W21,p->W22);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || coupling==0)
        continue;

        double *b = &buf[n*nv];

        if(owns(p,B.x(0),B.x(1),B.x(2)))
        {
            b[0] = p->ccipol1(a->u,B.x(0),B.x(1),B.x(2));
            b[1] = p->j_dir==1 ? p->ccipol2(a->v,B.x(0),B.x(1),B.x(2)) : 0.0;
            b[2] = p->ccipol3(a->w,B.x(0),B.x(1),B.x(2));
            b[3] = p->ccipol4a(a->ro,B.x(0),B.x(1),B.x(2));
            b[4] = p->ccipol4a(a->visc,B.x(0),B.x(1),B.x(2));
            b[5] = 1.0 - p->ccipol4a(*ALPHA,B.x(0),B.x(1),B.x(2));
            b[6] = 1.0;
        }

        if(B.mode==0)
        {
            const dem_shape &S = core.shapes[B.shape];
            for(size_t q=0; q<S.qp.size(); ++q)
            {
                dem_vec r = B.R*S.qp[q];
                dem_vec xq = B.x + r;
                if(!owns(p,xq(0),xq(1),xq(2)))
                continue;

                double ro = p->ccipol4a(a->ro,xq(0),xq(1),xq(2));
                dem_vec F = -S.qw[q]*ro*g;
                dem_vec T = r.cross(F);
                b[7]+=F(0); b[8]+=F(1); b[9]+=F(2);
                b[10]+=T(0); b[11]+=T(1); b[12]+=T(2);
                b[13]+=S.qw[q];
            }
        }
    }

    reduce_owner(pgc,buf,nv,false,false);

    for(int n=0; n<nb; ++n)
    {
        const double *b = &buf[n*nv];
        int cnt = int(b[6]+0.5);

        core.bodies[n].cpl.ufl_old = core.bodies[n].cpl.ufl;
        bool wasvalid = core.bodies[n].cpl.fluidcount>0;
        core.bodies[n].cpl.fluidcount = cnt;

        if(cnt>0)
        {
            core.bodies[n].cpl.ufl = dem_vec(b[0],b[1],b[2])/double(cnt);
            core.bodies[n].cpl.rhof = b[3]/double(cnt);
            core.bodies[n].cpl.nuf = b[4]/double(cnt);
            core.bodies[n].cpl.epsf = std::max(0.0,std::min(1.0,b[5]/double(cnt)));
        }
        core.bodies[n].cpl.ufl_valid = wasvalid && cnt>0;

        core.bodies[n].cpl.Fb = dem_vec(b[7],b[8],b[9]);
        core.bodies[n].cpl.Tb = dem_vec(b[10],b[11],b[12]);
        core.bodies[n].cpl.vsub = b[13];
    }
}

// ---------------------------------------------------------------------------------------------
// unresolved particles: reaction force and solid fraction spread to the grid
// ---------------------------------------------------------------------------------------------

void dem_f::feedback_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    ULOOP
    (*Sx)(i,j,k) = 0.0;
    VLOOP
    (*Sy)(i,j,k) = 0.0;
    WLOOP
    (*Sz)(i,j,k) = 0.0;
    LOOP
    (*ALPHA)(i,j,k) = 0.0;

    if(coupling==0 || coupling==2)
    return;

    int i0,i1,j0,j1,k0,k1;
    vector<double> sw(4*nb,0.0);

    // pass 1: kernel sums
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || B.mode!=0)
        continue;

        double R = std::max(core.shapes[B.shape].deq,kernel_cells*dxs);
        cellrange(p,n,R,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            if(p->flag1[IJK]>0)
            sw[4*n+0] += kernel(relpos(p,p->pos1_x(),p->pos1_y(),p->pos1_z(),B.x).norm(),R)*p->DXP[IP]*p->DYN[JP]*p->DZN[KP];
            if(p->flag2[IJK]>0)
            sw[4*n+1] += kernel(relpos(p,p->pos2_x(),p->pos2_y(),p->pos2_z(),B.x).norm(),R)*p->DXN[IP]*p->DYP[JP]*p->DZN[KP];
            if(p->flag3[IJK]>0)
            sw[4*n+2] += kernel(relpos(p,p->pos3_x(),p->pos3_y(),p->pos3_z(),B.x).norm(),R)*p->DXN[IP]*p->DYN[JP]*p->DZP[KP];
            if(p->flag4[IJK]>0)
            sw[4*n+3] += kernel(relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x).norm(),R)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        }
    }

    // kernel sums at the owners, then back to the ghosts
    reduce_owner(pgc,sw,4,true,false);

    // reaction force (drag and added mass) at the owners, sent to the ghosts
    vector<double> ffp(3*nb,0.0);
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(B.ghost || !B.active || B.fixed || B.mode!=0)
        continue;
        dem_vec ap = (B.v - B.cpl.vprev)/p->dt;
        B.cpl.Ffp = -(B.K*(B.uf - B.v) + B.madd*(B.af - ap));
        for(int q=0; q<3; ++q)
        ffp[3*n+q] = B.cpl.Ffp(q);
    }
    owner_to_ghosts(pgc,ffp,3);
    for(int n=0; n<nb; ++n)
    if(core.bodies[n].ghost)
    core.bodies[n].cpl.Ffp = dem_vec(ffp[3*n],ffp[3*n+1],ffp[3*n+2]);

    // pass 2: spread the reaction force and the solid volume on every copy of the particle
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || B.mode!=0)
        continue;

        const dem_shape &S = core.shapes[B.shape];

        double R = std::max(S.deq,kernel_cells*dxs);
        cellrange(p,n,R,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            double wk;

            if(p->flag1[IJK]>0 && sw[4*n+0]>0.0)
            {
                wk = kernel(relpos(p,p->pos1_x(),p->pos1_y(),p->pos1_z(),B.x).norm(),R);
                double ro = 0.5*(a->ro(i,j,k)+a->ro(i+1,j,k));
                (*Sx)(i,j,k) += core.bodies[n].cpl.Ffp(0)*wk/(ro*sw[4*n+0]);
            }

            if(p->flag2[IJK]>0 && sw[4*n+1]>0.0 && p->j_dir==1)
            {
                wk = kernel(relpos(p,p->pos2_x(),p->pos2_y(),p->pos2_z(),B.x).norm(),R);
                double ro = 0.5*(a->ro(i,j,k)+a->ro(i,j+1,k));
                (*Sy)(i,j,k) += core.bodies[n].cpl.Ffp(1)*wk/(ro*sw[4*n+1]);
            }

            if(p->flag3[IJK]>0 && sw[4*n+2]>0.0)
            {
                wk = kernel(relpos(p,p->pos3_x(),p->pos3_y(),p->pos3_z(),B.x).norm(),R);
                double ro = 0.5*(a->ro(i,j,k)+a->ro(i,j,k+1));
                (*Sz)(i,j,k) += core.bodies[n].cpl.Ffp(2)*wk/(ro*sw[4*n+2]);
            }

            if(p->flag4[IJK]>0 && sw[4*n+3]>0.0)
            {
                wk = kernel(relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x).norm(),R);
                (*ALPHA)(i,j,k) += S.volume*wk/sw[4*n+3];
            }
        }
    }

    LOOP
    (*ALPHA)(i,j,k) = std::min((*ALPHA)(i,j,k),0.8);

    pgc->start4a(p,*ALPHA,1);
}

// ---------------------------------------------------------------------------------------------
// forcing inside the RK stages of the CFD momentum step
// ---------------------------------------------------------------------------------------------

void dem_f::forcing_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, double alpha, field &u, field &v, field &w, bool finalize)
{
    if(!initialized || coupling==0)
    return;

    double starttime = pgc->timer();

    // unresolved: momentum source
    if(coupling==1 || coupling==3)
    {
        ULOOP
        u(i,j,k) += alpha*p->dt*(*Sx)(i,j,k);

        if(p->j_dir==1)
        VLOOP
        v(i,j,k) += alpha*p->dt*(*Sy)(i,j,k);

        WLOOP
        w(i,j,k) += alpha*p->dt*(*Sz)(i,j,k);
    }

    // resolved: direct forcing of the rigid body velocity inside the particle
    int i0,i1,j0,j1,k0,k1;
    double eps = hs_factor*dxs;
    int st = std::min(std::max(iter,0),2);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.mode!=1 || coupling==1)
        continue;

        const dem_shape &S = core.shapes[B.shape];
        double rlim = S.rbound + eps;
        cellrange(p,n,rlim+dxs,i0,i1,j0,j1,k0,k1);

        auto forcepoint = [&](double xf, double yf, double zf, int comp, double uval, double &fval, double &Hval)
        {
            dem_vec r = relpos(p,xf,yf,zf,B.x);
            fval = 0.0;
            Hval = 0.0;
            if(r.squaredNorm()>rlim*rlim)
            return;
            Hval = heaviside(S.sdf(B.R.transpose()*r),eps);
            if(Hval<=0.0)
            return;
            dem_vec ub = B.v + B.w.cross(r);
            fval = Hval*(ub(comp) - uval)/(alpha*p->dt);
        };

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            double f,H;

            if(p->flag1[IJK]>0)
            {
                double xf=p->pos1_x(), yf=p->pos1_y(), zf=p->pos1_z();
                forcepoint(xf,yf,zf,0,u(i,j,k),f,H);
                if(H>0.0)
                {
                    u(i,j,k) += alpha*p->dt*f;
                    double dF = -0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*f*p->DXP[IP]*p->DYN[JP]*p->DZN[KP];
                    core.bodies[n].cpl.Fs[st](0) += dF;
                    core.bodies[n].cpl.Ts[st] += relpos(p,xf,yf,zf,B.x).cross(dem_vec(dF,0.0,0.0));
                }
            }

            if(p->flag2[IJK]>0 && p->j_dir==1)
            {
                double xf=p->pos2_x(), yf=p->pos2_y(), zf=p->pos2_z();
                forcepoint(xf,yf,zf,1,v(i,j,k),f,H);
                if(H>0.0)
                {
                    v(i,j,k) += alpha*p->dt*f;
                    double dF = -0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*f*p->DXN[IP]*p->DYP[JP]*p->DZN[KP];
                    core.bodies[n].cpl.Fs[st](1) += dF;
                    core.bodies[n].cpl.Ts[st] += relpos(p,xf,yf,zf,B.x).cross(dem_vec(0.0,dF,0.0));
                }
            }

            if(p->flag3[IJK]>0)
            {
                double xf=p->pos3_x(), yf=p->pos3_y(), zf=p->pos3_z();
                forcepoint(xf,yf,zf,2,w(i,j,k),f,H);
                if(H>0.0)
                {
                    w(i,j,k) += alpha*p->dt*f;
                    double dF = -0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*f*p->DXN[IP]*p->DYN[JP]*p->DZP[KP];
                    core.bodies[n].cpl.Fs[st](2) += dF;
                    core.bodies[n].cpl.Ts[st] += relpos(p,xf,yf,zf,B.x).cross(dem_vec(0.0,0.0,dF));
                }
            }

            if(finalize && p->flag4[IJK]>0)
            {
                dem_vec r = relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x);
                if(r.squaredNorm()<=rlim*rlim)
                {
                    double Hc = heaviside(S.sdf(B.R.transpose()*r),eps);
                    core.bodies[n].cpl.mfl += a->ro(i,j,k)*Hc*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                    core.bodies[n].cpl.hvol += Hc*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                }
            }
        }
    }

    // global condition: combine_stages communicates
    if(finalize && (coupling==2 || coupling==3))
    combine_stages(p,pgc,iter,alpha);

    pgc->start1(p,u,10);
    pgc->start2(p,v,11);
    pgc->start3(p,w,12);

    p->fbtime += pgc->timer()-starttime;
}

// ---------------------------------------------------------------------------------------------
// contacts with the REEF3D topo and solid level sets
// ---------------------------------------------------------------------------------------------

double dem_f::wallphi_cfd(lexer *p, fdm *a, double x, double y, double z)
{
    double phi = 1.0e20;

    if(p->toporead>0 || p->S10>0)
    phi = std::min(phi,p->ccipol4a(a->topo,x,y,z));

    if(p->solidread==1)
    phi = std::min(phi,p->ccipol4a(a->solid,x,y,z));

    return phi;
}

void dem_f::walls_cfd(lexer *p, fdm *a, ghostcell *pgc, double margin, vector<dem_contact> &cts)
{
    // wall distance at the particle centres for a quick rejection
    vector<double> phic(nb,1.0e20);
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(B.active && !B.fixed && owns(p,B.x(0),B.x(1),B.x(2)))
        phic[n] = wallphi_cfd(p,a,B.x(0),B.x(1),B.x(2));
    }
    reduce_owner(pgc,phic,1,true,true);

    double h = 0.5*dxs;
    vector<double> loc;

    auto normal = [&](const dem_vec &x)
    {
        dem_vec g;
        g(0) = wallphi_cfd(p,a,x(0)+h,x(1),x(2)) - wallphi_cfd(p,a,x(0)-h,x(1),x(2));
        g(1) = p->j_dir==1 ? wallphi_cfd(p,a,x(0),x(1)+h,x(2)) - wallphi_cfd(p,a,x(0),x(1)-h,x(2)) : 0.0;
        g(2) = wallphi_cfd(p,a,x(0),x(1),x(2)+h) - wallphi_cfd(p,a,x(0),x(1),x(2)-h);
        double nn = g.norm();
        return nn>1.0e-14 ? dem_vec(g/nn) : dem_vec(dem_vec::UnitZ());
    };

    auto push = [&](int n, int feature, const dem_vec &x, const dem_vec &nrm, double gap)
    {
        loc.insert(loc.end(),{double(n),double(feature),x(0),x(1),x(2),nrm(0),nrm(1),nrm(2),gap});
    };

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        const dem_shape &S = core.shapes[B.shape];

        if(!B.active || B.fixed || phic[n]>1.0e19 || phic[n]-S.rbound>margin)
        continue;

        if(S.type==DEM_SPHERE)
        {
            if(owns(p,B.x(0),B.x(1),B.x(2)))
            {
                double gap = phic[n] - S.dim(0);
                if(gap<margin)
                {
                    dem_vec nrm = normal(B.x);
                    push(n,0,B.x - nrm*(S.dim(0)+0.5*gap),nrm,gap);
                }
            }
            continue;
        }

        for(size_t q=0; q<S.nodes.size(); ++q)
        {
            dem_vec xw = B.x + B.R*S.nodes[q];
            if(!owns(p,xw(0),xw(1),xw(2)))
            continue;

            double phi = wallphi_cfd(p,a,xw(0),xw(1),xw(2));
            if(phi<margin)
            {
                dem_vec nrm = normal(xw);
                push(n,q,xw - 0.5*phi*nrm,nrm,phi);
            }
        }
    }

    route_walls(p,pgc,loc,cts);
}

// ---------------------------------------------------------------------------------------------
// resolved particles: fluid momentum and angular momentum inside the smoothed particle indicator
// ---------------------------------------------------------------------------------------------

void dem_f::internal_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(coupling==0 || coupling==1)
    return;

    int i0,i1,j0,j1,k0,k1;
    double eps = hs_factor*dxs;
    vector<double> buf(6*nb,0.0);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.mode!=1 || B.fixed)
        continue;

        const dem_shape &S = core.shapes[B.shape];
        double rlim = S.rbound + eps;
        cellrange(p,n,rlim+dxs,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            if(p->flag4[IJK]<=0)
            continue;

            dem_vec r = relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x);
            if(r.squaredNorm()>rlim*rlim)
            continue;

            double H = heaviside(S.sdf(B.R.transpose()*r),eps);
            if(H<=0.0)
            continue;

            dem_vec uc(0.5*(a->u(i,j,k)+a->u(i-1,j,k)), 0.5*(a->v(i,j,k)+a->v(i,j-1,k)), 0.5*(a->w(i,j,k)+a->w(i,j,k-1)));
            if(p->j_dir==0)
            uc(1) = 0.0;

            dem_vec mom = a->ro(i,j,k)*H*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*uc;
            dem_vec ang = r.cross(mom);
            for(int q=0; q<3; ++q)
            {
                buf[6*n+q] += mom(q);
                buf[6*n+3+q] += ang(q);
            }
        }
    }

    reduce_owner(pgc,buf,6,false,false);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(B.ghost || !B.active || B.mode!=1 || B.fixed)
        continue;

        core.bodies[n].cpl.Ifl_old = core.bodies[n].cpl.Ifl;
        core.bodies[n].cpl.Lfl_old = core.bodies[n].cpl.Lfl;
        core.bodies[n].cpl.Ifl = dem_vec(buf[6*n],buf[6*n+1],buf[6*n+2]);
        core.bodies[n].cpl.Lfl = dem_vec(buf[6*n+3],buf[6*n+4],buf[6*n+5]);

        // first step: no history
        if(!core.bodies[n].cpl.Ifl_valid)
        {
            core.bodies[n].cpl.Ifl_old = core.bodies[n].cpl.Ifl;
            core.bodies[n].cpl.Lfl_old = core.bodies[n].cpl.Lfl;
            core.bodies[n].cpl.Ifl_valid = true;
        }
    }
}
