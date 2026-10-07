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

void sixdof_obj_nhflow::update_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // porous floating body: Darcy-Forchheimer resistance instead of rigid direct forcing
    if(p->X16==1)
    {
    update_forcing_nhflow_porous(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter);
    return;
    }
    
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
        
        uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
        vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
        wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
    }
    
    if(p->X15==1)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
        vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
        wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
    // Normal vectors calculation 
		nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
            - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        
        if(p->j_dir==0)
        ny = 0.0;
        
        if(p->j_dir==1)
		ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
            - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            
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
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        const double sb = tb*tb*(3.0 - 2.0*tb);
        
        FX[IJK] += (sb + (1.0-sb)*fabs(nx))*H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += (sb + (1.0-sb)*fabs(ny))*H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += (sb + (1.0-sb)*fabs(nz))*H*(wf - W[IJK])/(alpha[iter]*p->dt);
    
    }
    
    
    if(p->X15==2)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
        vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
        wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
        
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
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
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
   
    pgc->start5V(p,d->FHB,50);
}
    
double sixdof_obj_nhflow::Hpsi_nhflow(lexer *p, fdm_nhf *d)
{
    // half width of the smoothed Heaviside of the direct forcing (as in Hsolidface_nhflow)
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
