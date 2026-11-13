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

#include"nhflow_sigma.h"
#include"lexer.h"
#include"nhflow_curv.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

#define WLVL (fabs(d->WL(i,j))>0.00005?d->WL(i,j):1.0e20)
#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)

void nhflow_sigma::omega_update(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *U, double *V, double *W)
{ 
    FLOOP
    d->omegaF[FIJK] = 0.0;
        
        
    if(p->geo_curv==1)
    LOOP
    d->omegaF[FIJKp1] =   d->omegaF[FIJK]
                        - p->DZN[KP]*(d->detadt(i,j) + nhflow_curv::divh(p,d->FEx,d->FEy,i,j,k));
    
    if(p->geo_curv==0)
    LOOP
    {
    d->omegaF[FIJKp1] =   d->omegaF[FIJK]
                            
                        - p->DZN[KP]*(d->detadt(i,j) 
                            
                        + (d->FEx[IJK] - d->FEx[Im1JK])/p->DXN[IP]  + (d->FEy[IJK] - d->FEy[IJm1K])/p->DYN[JP]*p->y_dir);
    }
    
    // VRANS (B 200 1, B 221 1): the sigma faces move with the free surface (z_t = sigma D_t) through the
    // porous structure, which is fixed in space. The water volume a moving face sweeps is n_f sigma D_t,
    // not sigma D_t (Reynolds transport of the fluid fraction, d/dt int_V(t) n dV = int n v_g.N dA):
    //      omega_f = Q_f - n_f sigma_f D_t,   Q_f the Darcy flux through the face from the layer continuity.
    // The loops above give Q_f - sigma_f D_t (storage 1 in every layer); with D_t = -sum(div)/n_s
    // (continuity n_s deta/dt + div Q = 0) the flux through the surface was sum(div)(1/n_s - 1) instead of 0,
    // and every porous layer below a moving surface got a spurious vertical flux (1-n_f) sigma_f D_t
    // (0.5 m/s at the run-down over an armour layer). Face porosity n_f: mean of the two cells, the
    // surface porosity n_s (d->PORS) at the free surface. Open water (n = 1) unchanged.
    if(p->B200==1 && p->B221==1)
    LOOP
    {
        const double nf = (k==p->knoz-1) ? d->PORS(i,j) : 0.5*(d->POR[IJK] + d->POR[IJKp1]);
        
        d->omegaF[FIJKp1] += (1.0 - nf)*p->ZN[KP1]*d->detadt(i,j);
    }
    
      
    GC4LOOP
    if(p->gcb4[n][3]==6 && p->gcb4[n][4]==3)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
        k+=1;
        d->omegaF[FIJK] =  0.0;
        d->omegaF[FIJKp1] =  0.0;
        d->omegaF[FIJKp2] =  0.0;
        d->omegaF[FIJKp3] =  0.0;
    }
    
    GC4LOOP
    if(p->gcb4[n][3]==5 && p->gcb4[n][4]==21)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
        
        d->omegaF[FIJK] =  0.0;
        d->omegaF[FIJKm1] =  0.0;
        d->omegaF[FIJKm2] =  0.0;
        d->omegaF[FIJKm3] =  0.0;
    }
    

    FLOOP
    if(p->wet[IJ]==0)
    d->omegaF[FIJK] = 0.0;
    
    pgc->start7S(p,d->omegaF,17);
    
    // ZSN
    FLOOP
    p->ZSN[FIJK] = p->ZN[KP]*WL(i,j) + d->bed(i,j);
    
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*WL(i,j) + d->bed(i,j);
    
    pgc->start7S(p,p->ZSN,1);
    pgc->start5V(p,p->ZSP,1);
    
    SLICELOOP4
    {
        k=0;

            p->ZSN[FIJKm1] = p->ZN[KM1]*WL(i,j) + d->bed(i,j);
            p->ZSN[FIJKm2] = p->ZN[KM2]*WL(i,j) + d->bed(i,j);
            p->ZSN[FIJKm3] = p->ZN[KM3]*WL(i,j) + d->bed(i,j);
        
        k=p->knoz;

            p->ZSN[FIJKp1] = p->ZN[KP1]*WL(i,j) + d->bed(i,j);
            p->ZSN[FIJKp2] = p->ZN[KP2]*WL(i,j) + d->bed(i,j);
            p->ZSN[FIJKp3] = p->ZN[KP3]*WL(i,j) + d->bed(i,j);
    }
    
    
    SLICELOOP4
    {
        k=0;

            p->ZSP[IJKm1] = p->ZP[KM1]*WL(i,j) + d->bed(i,j);
            p->ZSP[IJKm2] = p->ZP[KM2]*WL(i,j) + d->bed(i,j);
            p->ZSP[IJKm3] = p->ZP[KM3]*WL(i,j) + d->bed(i,j);
        
        k=p->knoz-1;

            p->ZSP[IJKp1] = p->ZP[KP1]*WL(i,j) + d->bed(i,j);
            p->ZSP[IJKp2] = p->ZP[KP2]*WL(i,j) + d->bed(i,j);
            p->ZSP[IJKp3] = p->ZP[KP3]*WL(i,j) + d->bed(i,j);
    }
    
}




