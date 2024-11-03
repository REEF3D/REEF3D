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
                             double alpha, double *U, double *V, double *W, double *FX, double *FY, double *FZ)
{/*
    // Calculate forcing fields
    double H,Ht, uf, vf, wf;
	double nx, ny, nz,norm ;
	double psi, phival_fb;
    double dirac;
    
    H=Ht=0.0;
  
    ULOOP
    {
        uf = u_fb(0) + u_fb(4)*(p->pos1_z() - c_(2)) - u_fb(5)*(p->pos1_y() - c_(1));
        
        H = Hsolidface(p,a,1,0,0);
        
        fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);   
        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    
    VLOOP
    {
        vf = u_fb(1) + u_fb(5)*(p->pos2_x() - c_(0)) - u_fb(3)*(p->pos2_z() - c_(2));
        
        H = Hsolidface(p,a,0,1,0);

        fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);
        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H, 1.0); 
    }
    
    WLOOP
    {
        wf = u_fb(2) + u_fb(3)*(p->pos3_y() - c_(1)) - u_fb(4)*(p->pos3_x() - c_(0));

        H = Hsolidface(p,a,0,0,1);
        
        fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);

    }
    
    LOOP
    {
        H = Hsolidface(p,a,0,0,0);
        a->fbh4(i,j,k) = min(a->fbh4(i,j,k) + H, 1.0); 
    }
 */   
}
    
    
double sixdof_obj::Hsolidface_nhflow(lexer *p, fdm_nhf *d, int aa, int bb, int cc)
{
    double psi, H, phival_fb,dirac;
    
    if (p->j_dir==0)
    psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]/p->sigz[IJ]);
	
    if (p->j_dir==1)
    psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]/p->sigz[IJ]);


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
    
    
    
