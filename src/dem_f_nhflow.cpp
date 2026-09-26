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
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<mpi.h>

// ---------------------------------------------------------------------------------------------
// fluid data at the particles: velocity and voidage at the centroid, submerged volume and
// hydrostatic buoyancy from the volume quadrature against the NHFLOW free surface
// ---------------------------------------------------------------------------------------------

void dem_f::fluid_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const int nv = 14;
    vector<double> buf(nb*nv,0.0);
    dem_vec g(p->W20,p->W21,p->W22);

    // resolved particles force the fluid inside them, including the free surface in their footprint;
    // their hydrostatic buoyancy therefore uses the ambient water level, sampled on a ring around them
    vector<double> ring(2*nb,0.0);
    const int nring = 24;
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || coupling==0 || core.bodies[n].cpl.basemode!=1)
        continue;

        double rr = core.shapes[B.shape].rbound + 2.0*dxs;
        for(int m=0; m<nring; ++m)
        {
            double phi = 2.0*3.14159265358979323846*m/nring;
            double xr = B.x(0) + rr*cos(phi);
            double yr = p->j_dir==1 ? B.x(1) + rr*sin(phi) : B.x(1);
            if(p->j_dir==0 && m%(nring/2)!=0)
            continue;
            if(!owns(p,xr,yr,B.x(2)))
            continue;
            ring[2*n]   += p->ccslipol4(d->WL,xr,yr) + p->ccslipol4(d->bed,xr,yr);
            ring[2*n+1] += 1.0;
        }
    }
    reduce_owner(pgc,ring,2,true,false);

    // NHFLOW has no free surface inside a forced particle: resolved forcing is used only while the
    // particle is fully submerged, surface-piercing particles are treated as unresolved
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(B.ghost || core.bodies[n].cpl.basemode!=1 || ring[2*n+1]<0.5)
        continue;

        double eta = ring[2*n]/ring[2*n+1];
        int mode = (B.x(2) + core.shapes[B.shape].rbound < eta - dxs) ? 1 : 0;
        if(mode!=B.mode)
        {
            B.mode = mode;
            core.bodies[n].cpl.Ifl_valid = false;
            if((B.tier==1 && p->mpirank==0) || (B.tier==0 && B.owner==p->mpirank))
            cout<<"DEM: particle "<<B.id<<(mode==1 ? " submerged, resolved coupling" : " surface-piercing, unresolved coupling")<<endl;
        }
    }

    // ghosts need the new mode before the forcing and the internal momentum of this step
    {
        vector<double> md(2*nb,0.0);
        for(int n=0; n<nb; ++n)
        {
            md[2*n]   = core.bodies[n].mode;
            md[2*n+1] = core.bodies[n].cpl.Ifl_valid ? 1.0 : 0.0;
        }
        owner_to_ghosts(pgc,md,2);
        for(int n=0; n<nb; ++n)
        if(core.bodies[n].ghost)
        {
            core.bodies[n].mode = int(llround(md[2*n]));
            core.bodies[n].cpl.Ifl_valid = md[2*n+1]>0.5;
        }
    }

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || coupling==0)
        continue;

        double *b = &buf[n*nv];

        if(owns(p,B.x(0),B.x(1),B.x(2)))
        {
            b[0] = p->ccipol4V(d->U,d->WL,d->bed,B.x(0),B.x(1),B.x(2));
            b[1] = p->j_dir==1 ? p->ccipol4V(d->V,d->WL,d->bed,B.x(0),B.x(1),B.x(2)) : 0.0;
            b[2] = p->ccipol4V(d->W,d->WL,d->bed,B.x(0),B.x(1),B.x(2));
            b[3] = p->W1;
            b[4] = p->W2;
            b[5] = 1.0 - p->ccipol4V(ALPHAV,d->WL,d->bed,B.x(0),B.x(1),B.x(2));
            b[6] = 1.0;
        }

        // large particles (resolved or surface-piercing): ambient water level, avoids self-interaction
        bool useRing = core.bodies[n].cpl.basemode==1 && ring[2*n+1]>0.5;
        double etaring = useRing ? ring[2*n]/ring[2*n+1] : 0.0;

        const dem_shape &S = core.shapes[B.shape];
        for(size_t q=0; q<S.qp.size(); ++q)
        {
            dem_vec r = B.R*S.qp[q];
            dem_vec xq = B.x + r;
            if(!owns(p,xq(0),xq(1),xq(2)))
            continue;

            double eta = useRing ? etaring : p->ccslipol4(d->WL,xq(0),xq(1)) + p->ccslipol4(d->bed,xq(0),xq(1));
            if(xq(2)>=eta)
            continue;

            dem_vec F = -S.qw[q]*p->W1*g;
            dem_vec T = r.cross(F);
            b[7]+=F(0); b[8]+=F(1); b[9]+=F(2);
            b[10]+=T(0); b[11]+=T(1); b[12]+=T(2);
            b[13]+=S.qw[q];
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
// unresolved particles: reaction force and solid fraction spread to the sigma grid
// ---------------------------------------------------------------------------------------------

void dem_f::feedback_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    {
        SX[IJK] = 0.0;
        SY[IJK] = 0.0;
        SZ[IJK] = 0.0;
        ALPHAV[IJK] = 0.0;
    }

    if(coupling==0 || coupling==2)
    return;

    int i0,i1,j0,j1,k0,k1;
    vector<double> sw(nb,0.0);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];

        // surface-piercing particles larger than the grid are coupled one-way (fluid to particle)
        if(!B.active || B.fixed || B.mode!=0 || core.bodies[n].cpl.basemode==1)
        continue;

        double R = std::max(core.shapes[B.shape].deq,kernel_cells*dxs);
        cellrange(p,n,R,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        if(p->flag4[IJK]>0 && p->wet[IJ]>0)
        sw[n] += kernel(relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x).norm(),R)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*d->WL(i,j);
    }

    reduce_owner(pgc,sw,1,true,false);

    // reaction force and submerged volume at the owners, sent to the ghosts
    vector<double> ffp(4*nb,0.0);
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(B.ghost || !B.active || B.fixed || B.mode!=0)
        continue;
        dem_vec ap = (B.v - B.cpl.vprev)/p->dt;
        B.cpl.Ffp = -(B.K*(B.uf - B.v) + B.madd*(B.af - ap));
        for(int q=0; q<3; ++q)
        ffp[4*n+q] = B.cpl.Ffp(q);
        ffp[4*n+3] = B.cpl.vsub;
    }
    owner_to_ghosts(pgc,ffp,4);
    for(int n=0; n<nb; ++n)
    if(core.bodies[n].ghost)
    {
        core.bodies[n].cpl.Ffp = dem_vec(ffp[4*n],ffp[4*n+1],ffp[4*n+2]);
        core.bodies[n].cpl.vsub = ffp[4*n+3];
    }

    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || B.mode!=0 || core.bodies[n].cpl.basemode==1 || sw[n]<=0.0)
        continue;

        const dem_shape &S = core.shapes[B.shape];

        double R = std::max(S.deq,kernel_cells*dxs);
        cellrange(p,n,R,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        if(p->flag4[IJK]>0 && p->wet[IJ]>0)
        {
            double wk = kernel(relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x).norm(),R);
            SX[IJK] += core.bodies[n].cpl.Ffp(0)*wk/(p->W1*sw[n]);
            SY[IJK] += p->j_dir==1 ? core.bodies[n].cpl.Ffp(1)*wk/(p->W1*sw[n]) : 0.0;
            SZ[IJK] += core.bodies[n].cpl.Ffp(2)*wk/(p->W1*sw[n]);
            ALPHAV[IJK] += core.bodies[n].cpl.vsub*wk/sw[n];
        }
    }

    LOOP
    ALPHAV[IJK] = std::min(ALPHAV[IJK],0.8);

    pgc->start4V(p,ALPHAV,1);
}

// ---------------------------------------------------------------------------------------------
// forcing inside the RK stages of the NHFLOW momentum step
// ---------------------------------------------------------------------------------------------

void dem_f::forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double alpha,
                           double *UH, double *VH, double *WH, slice &WL, bool finalize, bool reforce)
{
    if(!initialized || coupling==0)
    return;

    double starttime = pgc->timer();

    // unresolved: momentum source (once per stage, not in the re-forcing after the projection)
    if((coupling==1 || coupling==3) && !reforce)
    LOOP
    {
        d->U[IJK] += alpha*p->dt*SX[IJK];
        UH[IJK]   += alpha*p->dt*SX[IJK]*WL(i,j);

        d->V[IJK] += alpha*p->dt*SY[IJK];
        VH[IJK]   += alpha*p->dt*SY[IJK]*WL(i,j);

        d->W[IJK] += alpha*p->dt*SZ[IJK];
        WH[IJK]   += alpha*p->dt*SZ[IJK]*WL(i,j);
    }

    // resolved: direct forcing
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

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            if(p->flag4[IJK]<=0 || p->wet[IJ]==0)
            continue;

            dem_vec r = relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x);
            if(r.squaredNorm()>rlim*rlim)
            continue;

            double H = heaviside(S.sdf(B.R.transpose()*r),eps);
            if(H<=0.0)
            continue;

            dem_vec ub = B.v + B.w.cross(r);
            dem_vec f;
            f(0) = H*(ub(0) - d->U[IJK])/(alpha*p->dt);
            f(1) = p->j_dir==1 ? H*(ub(1) - d->V[IJK])/(alpha*p->dt) : 0.0;
            f(2) = H*(ub(2) - d->W[IJK])/(alpha*p->dt);

            d->U[IJK] += alpha*p->dt*f(0);
            d->V[IJK] += alpha*p->dt*f(1);
            d->W[IJK] += alpha*p->dt*f(2);
            UH[IJK] += alpha*p->dt*f(0)*WL(i,j);
            VH[IJK] += alpha*p->dt*f(1)*WL(i,j);
            WH[IJK] += alpha*p->dt*f(2)*WL(i,j);

            // forcing before the projection and re-forcing after it (the latter carries the
            // non-hydrostatic pressure impulse) both belong to the hydrodynamic force
            {
                double dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*WL(i,j);
                dem_vec dF = -p->W1*f*dV;
                core.bodies[n].cpl.Fs[st] += dF;
                core.bodies[n].cpl.Ts[st] += r.cross(dF);
                if(finalize && !reforce)
                {
                    core.bodies[n].cpl.mfl += p->W1*H*dV;
                    core.bodies[n].cpl.hvol += H*dV;
                }
            }
        }
    }

    // global condition: combine_stages communicates
    if(finalize && reforce && (coupling==2 || coupling==3))
    combine_stages(p,pgc,iter,alpha);

    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
    pgc->start4V(p,UH,14);
    pgc->start4V(p,VH,15);
    pgc->start4V(p,WH,16);

    p->dftime += pgc->timer()-starttime;
}

// ---------------------------------------------------------------------------------------------
// contacts with the NHFLOW bed and solid level set
// ---------------------------------------------------------------------------------------------

double dem_f::wallphi_nhflow(lexer *p, fdm_nhf *d, double x, double y, double z)
{
    // bed as a level set: vertical distance scaled by the slope
    double h = 0.5*dxs;
    double bed = p->ccslipol4(d->bed,x,y);
    double bx = (p->ccslipol4(d->bed,x+h,y) - p->ccslipol4(d->bed,x-h,y))/(2.0*h);
    double by = p->j_dir==1 ? (p->ccslipol4(d->bed,x,y+h) - p->ccslipol4(d->bed,x,y-h))/(2.0*h) : 0.0;
    double phi = (z-bed)/sqrt(1.0 + bx*bx + by*by);

    if(p->solidread==1)
    phi = std::min(phi,p->ccipol4V(d->SOLID,d->WL,d->bed,x,y,z));

    return phi;
}

void dem_f::walls_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double margin, vector<dem_contact> &cts)
{
    vector<double> phic(nb,1.0e20);
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(B.active && !B.fixed && owns(p,B.x(0),B.x(1),B.x(2)))
        phic[n] = wallphi_nhflow(p,d,B.x(0),B.x(1),B.x(2));
    }
    reduce_owner(pgc,phic,1,true,true);

    double h = 0.5*dxs;
    vector<double> loc;

    auto normal = [&](const dem_vec &x)
    {
        dem_vec g;
        g(0) = wallphi_nhflow(p,d,x(0)+h,x(1),x(2)) - wallphi_nhflow(p,d,x(0)-h,x(1),x(2));
        g(1) = p->j_dir==1 ? wallphi_nhflow(p,d,x(0),x(1)+h,x(2)) - wallphi_nhflow(p,d,x(0),x(1)-h,x(2)) : 0.0;
        g(2) = wallphi_nhflow(p,d,x(0),x(1),x(2)+h) - wallphi_nhflow(p,d,x(0),x(1),x(2)-h);
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

            double phi = wallphi_nhflow(p,d,xw(0),xw(1),xw(2));
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

void dem_f::internal_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
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
            if(p->flag4[IJK]<=0 || p->wet[IJ]==0)
            continue;

            dem_vec r = relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x);
            if(r.squaredNorm()>rlim*rlim)
            continue;

            double H = heaviside(S.sdf(B.R.transpose()*r),eps);
            if(H<=0.0)
            continue;

            dem_vec uc(d->U[IJK], p->j_dir==1 ? d->V[IJK] : 0.0, d->W[IJK]);
            dem_vec mom = p->W1*H*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*d->WL(i,j)*uc;
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

        if(!core.bodies[n].cpl.Ifl_valid)
        {
            core.bodies[n].cpl.Ifl_old = core.bodies[n].cpl.Ifl;
            core.bodies[n].cpl.Lfl_old = core.bodies[n].cpl.Lfl;
            core.bodies[n].cpl.Ifl_valid = true;
        }
    }
}
