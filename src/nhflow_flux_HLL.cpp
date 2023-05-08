/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_flux_HLL.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_flux_HLL::nhflow_flux_HLL(lexer* p, patchBC_interface *ppBC) : nhflow_flux_reconstruct(p,ppBC)
{
    pBC = ppBC;
    
    p->Darray(Fs,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fn,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fe,p->imax*p->jmax*(p->kmax+2));
    p->Darray(Fw,p->imax*p->jmax*(p->kmax+2));
}

nhflow_flux_HLL::~nhflow_flux_HLL()
{
}

void nhflow_flux_HLL::face_flux_2D(lexer* p, fdm_nhf*, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
{
}

void nhflow_flux_HLL::face_flux_3D(lexer *p, ghostcell *pgc, fdm_nhf *d, slice & eta, double *U, double *V, double *Fx, double *Fy)
{
    double Ss,Sn,Se,Sw;
    double USx,USy;
    double Ds,Dn,De,Dw;
    double DSx,DSy;
    double denom;

         
    // reconstruct eta
    reconstruct_2D(p, pgc, d, eta, d->ETAs, d->ETAn, d->ETAe, d->ETAw);
    
    // reconstruct U and V
    reconstruct_3D(p, pgc, d, U, V, Fs, Fn, Fe, Fw);
    
    // HLL flux
    ULOOP
    {
    // water level       
    Ds = d->ETAs(i,j) + 0.5*(d->depth(i,j) + d->depth(i+1,j));
    Dn = d->ETAn(i,j) + 0.5*(d->depth(i,j) + d->depth(i+1,j));
    
    Ds = MAX(0.00005, Ds);
    Dn = MAX(0.00005, Dn);
    
    // Us
    USx = 0.5*(Fs[IJK]+Fn[IJK]) + sqrt(9.81*Ds) - sqrt(9.81*Dn);
    DSx = 0.5*(sqrt(9.81*Ds) + sqrt(9.81*Dn)) + 0.25*(Fs[IJK] - Fn[IJK]);
    
    // wave speed
    Ss = MIN(Fs[IJK] - sqrt(9.81*Ds), USx - DSx);
    Sn = MAX(Fn[IJK] + sqrt(9.81*Dn), USx + DSx);
    
        // final flux x-dir
        if(Ss>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(Sn<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        {
        denom = Sn-Ss;
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fx[IJK] = (Sn*Fs[IJK] - Ss*Fn[IJK] + Sn*Ss*(Dn - Ds))/denom;
        }
    }
    
    
    VLOOP
    {
    // water level       
    De = d->ETAe(i,j) + 0.5*(d->depth(i,j) + d->depth(i,j+1));
    Dw = d->ETAw(i,j) + 0.5*(d->depth(i,j) + d->depth(i,j+1));
    
    De = MAX(0.00005, De);
    Dw = MAX(0.00005, Dw);
    
    // Us
    USy = 0.5*(Fe[IJK]+Fw[IJK]) + sqrt(9.81*De) - sqrt(9.81*Dw);
    DSy = 0.5*(sqrt(9.81*De) + sqrt(9.81*Dw)) + 0.25*(Fe[IJK] - Fw[IJK]);
    
    // wave speed
    Se = MIN(Fe[IJK] - sqrt(9.81*De), USy - DSy);
    Sw = MAX(Fw[IJK] + sqrt(9.81*Dw), USy + DSy);

        
        // final flux y-dir
        if(Se>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(Sw<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        {
        denom = Sw-Se;
        denom = fabs(denom)>1.0e-10?denom:1.0e10;
        
        Fy[IJK] = (Sw*Fe[IJK] - Se*Fw[IJK] + Sw*Se*(Dw - De))/denom;
        }
    }
    
    pgc->start1V(p,Fx,10);
    pgc->start2V(p,Fy,10);
	
}


