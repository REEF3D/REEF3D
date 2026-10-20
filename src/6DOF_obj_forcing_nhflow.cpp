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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans_definitions.h"

void sixdof_obj_nhflow::update_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // porous floating body: Darcy-Forchheimer resistance instead of rigid direct forcing
    if(p->X16==1)
    {
    update_forcing_nhflow_porous(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter);
    return;
    }
    
    forcing_kernel_nhflow(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,false);
}

void sixdof_obj_nhflow::update_forcing_nhflow_wavemaker(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // wavemaker (X 10 4): the same direct forcing (X 15 0/1/2) as the floating body, with the
    // paddle velocity of the stage as the target
    forcing_kernel_nhflow(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,true);
}

void sixdof_obj_nhflow::target_velocity(lexer *p, bool wavemaker, double &uf, double &vf, double &wf)
{
    // target velocity of the direct forcing at the cell (i,j,k)
    if(wavemaker)
    {
        uf = wm_u(p,p->pos_x(),p->pos_z());
        vf = 0.0;
        wf = 0.0;
        return;
    }
    
    uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
    vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
    wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
}

void sixdof_obj_nhflow::forcing_kernel_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, int iter, bool wavemaker)
{
    // Calculate forcing fields
    double H, uf, vf, wf;
    double ef,efc;
    
    double du,dv,dw,udotn;
    double beta;   // 0.0 = free slip (tangential preserved), 1.0 = no slip
    
    if(p->X15==0)
    beta=1.0;
    
    if(p->X15>=1)
    beta=0.0;
    
    
    if(p->X15==0)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        target_velocity(p,wavemaker,uf,vf,wf);
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
    }
    
    // X 15 3: as X 15 1, but only the cells inside the body are forced. A body moving through the
    // grid at speed (towed or self-propelled ship) otherwise pushes the water in the forcing band ahead
    // of the bow up to the body velocity; that water passes under the hull and is carried along with it.
    // X 15 5: as X 15 1 with the sharp switch at FB = 0 (full forcing for FB <= 0, normal-weighted outside),
    // the forcing before the shell ramp was introduced (9ce616dab)
    // X 15 6 / 7: the forcing of X 15 1 / 3; the body faces are reflective walls in the HLL fluxes (nhflow_HLL.cpp)
    if(p->X15==1 || p->X15==3 || p->X15==5 || p->X15==6 || p->X15==7)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        target_velocity(p,wavemaker,uf,vf,wf);
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        if((p->X15==3 || p->X15==7) && d->FB[IJK]>=0.0)
        continue;
        
    // Normal vectors calculation 
        if(p->geo_curv==1)
        fbnormal_curv(p,d,nx,ny);
        else
        {
		nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
            - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        
        if(p->j_dir==0)
        ny = 0.0;
        
        if(p->j_dir==1)
		ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
            - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        }
            
		nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        
        // full forcing inside the body, normal-weighted forcing (tangential velocity kept) outside,
        // ramped over a thin shell |FB| < dsh around the surface. A sharp switch at FB = 0 jumps
        // for cell centres on or next to the surface (a flat side through a row of cell centres):
        // the smallest rotation moves the row partly inside and partly outside, and the jump in
        // the forcing kicks a symmetric hull into sway and yaw.
        // X 15 5: the sharp switch at FB = 0 (2D cases and grid-aligned hulls, where the kick does not
        // arise; reproduces the floating-body motions of the code before the ramp)
        double sb;
        
        if(p->X15==5)
        sb = d->FB[IJK]<=0.0 ? 1.0 : 0.0;
        
        else
        {
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        sb = tb*tb*(3.0 - 2.0*tb);
        }
        
        FX[IJK] += (sb + (1.0-sb)*fabs(nx))*H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += (sb + (1.0-sb)*fabs(ny))*H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += (sb + (1.0-sb)*fabs(nz))*H*(wf - W[IJK])/(alpha[iter]*p->dt);
    
    }
    
    
    if(p->X15==2)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        du = uf - U[IJK];
        dv = vf - V[IJK];
        dw = wf - W[IJK];
        
        if(d->FB[IJK] < 0.0)
        {
            // interior of the body: enforce the full rigid-body velocity
            FX[IJK] += H*du/(alpha[iter]*p->dt);
            FY[IJK] += H*dv/(alpha[iter]*p->dt);
            FZ[IJK] += H*dw/(alpha[iter]*p->dt);
        }
         
        else
        {
            // Normal vectors calculation 
            if(p->geo_curv==1)
            fbnormal_curv(p,d,nx,ny);
            else
            {
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            }
                
            nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

            norm = sqrt(nx*nx + ny*ny + nz*nz);
                    
            if(norm > 1.0e-10)
            {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            udotn = nx*du + ny*dv + nz*dw;
            
            FX[IJK] += H*(beta*du + (1.0-beta)*udotn*nx)/(alpha[iter]*p->dt);
            FY[IJK] += H*(beta*dv + (1.0-beta)*udotn*ny)/(alpha[iter]*p->dt);
            FZ[IJK] += H*(beta*dw + (1.0-beta)*udotn*nz)/(alpha[iter]*p->dt);
            }
            
            else
            {
            // degenerate gradient (flat level set): fall back to full forcing
            FX[IJK] += H*du/(alpha[iter]*p->dt);
            FY[IJK] += H*dv/(alpha[iter]*p->dt);
            FZ[IJK] += H*dw/(alpha[iter]*p->dt);
            }
        }
    
    }
    
    // X 15 4: normal projection outside the body (as X 15 2), ramped over the thin shell |FB| < dsh of X 15 1:
    //   F = H [sb du + (1 - sb)(du.n) n] / (alpha dt CPORNH)
    // sb = 1 inside the body (full rigid-body velocity), 0 outside (only the normal velocity is imposed,
    // the tangential velocity is kept). Unlike the component weighting |n_i| of X 15 1, the projection does
    // not depend on the orientation of the hull to the grid (rolled hull, bilges): the tangential velocity
    // outside is not damped. Unlike the sharp switch of X 15 2 at FB = 0, the forcing stays continuous when
    // the hull moves across cell centres. On grid-aligned faces X 15 1 and 4 are identical.
    // The forcing is divided by CPORNH as in nhflow_forcing::solid_forcing, so that the body velocity is
    // fully imposed also inside a porous structure (without porous media CPORNH = 1).
    if(p->X15==4)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        if(H<=0.0)
        continue;
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        du = uf - U[IJK];
        dv = vf - V[IJK];
        dw = wf - W[IJK];
        
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        double sb = tb*tb*(3.0 - 2.0*tb);
        
        nx=ny=nz=0.0;
        udotn=0.0;
        
        if(sb<1.0)
        {
            // Normal vectors calculation 
            if(p->geo_curv==1)
            fbnormal_curv(p,d,nx,ny);
            else
            {
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            }
                
            nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            if(norm > 1.0e-10)
            {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            udotn = nx*du + ny*dv + nz*dw;
            }
            
            // degenerate gradient (flat level set): fall back to full forcing
            else
            sb = 1.0;
        }
        
        const double fac = H/(alpha[iter]*p->dt*CPORNH);
        
        FX[IJK] += fac*(sb*du + (1.0-sb)*udotn*nx);
        FY[IJK] += fac*(sb*dv + (1.0-sb)*udotn*ny);
        FZ[IJK] += fac*(sb*dw + (1.0-sb)*udotn*nz);
    }
   
    pgc->start5V(p,d->FHB,50);
}
    
void sixdof_obj_nhflow::fbnormal_curv(lexer *p, fdm_nhf *d, double &gx, double &gy)
{
    // curvilinear grid (G 50): the horizontal part of -grad(FB) at the cell (i,j,k), the gradient along the sigma
    // level by the chain rule on the cell centres, then the sigma correction (sigx, sigy are physical)
    double fx,fy;
    p->geo_grad(i,j,d->FB[Ip1JK],d->FB[Im1JK],d->FB[IJp1K],d->FB[IJm1K],fx,fy);
    const double fs = (d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);

    gx = -fx - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*fs;
    gy = -fy - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*fs;
}

double sixdof_obj_nhflow::Hpsi_nhflow(lexer *p, fdm_nhf *d)
{
    // half width of the smoothed Heaviside of the direct forcing (as in Hsolidface_nhflow)
    if(p->geo_curv==1)
    return p->A526*p->geo_width(i,j);

    if(p->j_dir==0)
    return p->A526*p->DXN[IP];
    
    return p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);
}

double sixdof_obj_nhflow::Hsolidface_nhflow(lexer *p, fdm_nhf *d, int aa, int bb, int cc)
{
    double psi, H, phival_fb,dirac;
    
    if(p->j_dir==0)
    psi = p->A526*(1.0/1.0)*(p->DXN[IP] + 0.0*p->DZN[KP]*p->WL[IJ]);
    
    if(p->j_dir==1)
    psi = p->A526*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP] + 0.0*p->DZN[KP]*p->WL[IJ]);

    // curvilinear grid (G 50): the physical mean width of the cell
    if(p->geo_curv==1)
    psi = p->A526*p->geo_width(i,j);


    // Construct solid heaviside function
    phival_fb = d->FB[IJK];
    
    if(-phival_fb > psi)
    H = 1.0;

    if(-phival_fb < -psi)
    H = 0.0;

    if(fabs(phival_fb)<=psi)
    H = 0.5*(1.0 + (-phival_fb)/psi + (1.0/PI)*sin((PI*(-phival_fb))/psi));
    
    return H;
}
