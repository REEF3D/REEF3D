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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::update_forcing_nhflow_wavemaker(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // Calculate forcing fields
    double H, uf, vf, wf, pf, fac;
    double ef,efc;
    
    double x1,x2;
    
    if(p->X15==0)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        fac=1.0;
        
        if(p->X172==1)
        {
        x1 = (p->X172_xs + 0.5*(p->X172_xe-p->X172_xs));
        x2 = (p->X172_xs + 0.75*(p->X172_xe-p->X172_xs));

        fac = (p->XP[IP]-x1)*(1.0/(x2-x1));
        
        if(p->XP[IP] > x2)
        fac=1.0;
        
        if(p->XP[IP] < x1)
        fac=0.0;
        }
        
        //cout <<"fac: "<<fac<<endl;
        
        uf = fac*uwm[k]; 
        
        vf = 0.0;
        wf = 0.0;
        pf = 0.0;
        //cout <<"UF: "<<uwm[k]<<endl;
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
        
        //d->P[FIJK] += H*(pf - d->P[FIJK]);
    }
    
    if(p->X15==1)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        uf = uwm[k];
        vf = 0.0;
        wf = 0.0;
         
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
    // Normal vectors calculation 
		nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXN[IP] + p->DXN[IM1]);
		ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYN[JP] + p->DYN[JM1]);
		nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZN[KP]*WL(i,j) + p->DZN[KM1]*WL(i,j));

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        
        if(d->FB[IJK]<=0.0)
        {
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
        }

        if(d->FB[IJK]>0.0)
        {
        FX[IJK] += fabs(nx)*H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += fabs(ny)*H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += fabs(nz)*H*(wf - W[IJK])/(alpha[iter]*p->dt);
        }
    
    }
   
    pgc->start5V(p,d->FHB,50);
}
    

    
    
