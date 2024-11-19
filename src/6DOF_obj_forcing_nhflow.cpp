/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

void sixdof_obj::update_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // Calculate forcing fields
    double H, uf, vf, wf;
    double ef,efc;
    
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        d->FHB[IJK] = min(d->FHB[IJK] + H, 1.0); 
        
        uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
        vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
        wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
         
        
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
    }
    
    LOOP
    d->test[IJK] = FZ[IJK];
    
    k=p->knoz-1;
     
    SLICELOOP4
    {
    H = Hsolidface_nhflow(p,d,0,0,0);
    
    ef = d->bed(i,j) + d->depth(i,j);
    
    
    if(d->SOLID[IJK]<0.0)
    {
        ef = 0.0;
        efc = 0.0;
        
        if(d->SOLID[Im1JK]>0.0)   
        {
        ef += WL(i-1,j);
        efc+=1.0;
        }
        
        if(d->SOLID[Ip1JK]>0.0)    
        {
        ef += WL(i+1,j);
        efc+=1.0;
        }

        if(d->SOLID[IJm1K]>0.0) 
        {
        ef += WL(i,j-1);
        efc+=1.0;
        }
        
        if(d->SOLID[IJp1K]>0.0)    
        {
        ef += WL(i,j+1);
        efc+=1.0;
        }
        
    if(efc>0.1)
    ef = ef/efc;
    
    if(efc<0.1)
    ef = d->bed(i,j) + d->depth(i,j);
    }
    
    fe(i,j) += H*(ef - WL(i,j))/(alpha[iter]*p->dt);
    }

    pgc->start5V(p,d->FHB,1);
}
    
double sixdof_obj::Hsolidface_nhflow(lexer *p, fdm_nhf *d, int aa, int bb, int cc)
{
    double psi, H, phival_fb,dirac;
    
    if (p->j_dir==0)
    psi = p->X41*(1.0/1.0)*(p->DXN[IP]);// + 0.0*p->DZN[KP]*d->WL(i,j));
	
    if (p->j_dir==1)
    psi = p->X41*(1.0/2.0)*(p->DXN[IP]+p->DYN[JP]);// + 0.0*p->DZN[KP]*d->WL(i,j));


    // Construct solid heaviside function
    phival_fb = d->FB[IJK];
    
	
    if (-phival_fb > psi)
    H = 1.0;

    else if (-phival_fb < -psi)
    H = 0.0;

    else
    H = 0.5*(1.0 + -phival_fb/psi + (1.0/PI)*sin((PI*-phival_fb)/psi));


    return H;
}
    
    
    
