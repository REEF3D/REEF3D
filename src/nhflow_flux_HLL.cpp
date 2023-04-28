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
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_flux_HLL::nhflow_flux_HLL(lexer* p, patchBC_interface *ppBC) : nhflow_flux_reconstruct(p,ppBC),
                                                                    ETAs(p),ETAn(p),ETAe(p),ETAw(p)
{
    pBC = ppBC;
}

nhflow_flux_HLL::~nhflow_flux_HLL()
{
}

void nhflow_flux_HLL::face_flux_2D(lexer* p, fdm_nhf*, slice& f, slice &fs, slice &fn, slice &fe, slice &fw)
{
}

void nhflow_flux_HLL::face_flux_3D(lexer *p, fdm_nhf *d, slice & eta, double *F, double *Fs, double *Fn, double *Fe, double *Fw,
                                                                   double *Fx, double *Fy)
{
    double Ss,Sn,Se,Sw;
    double USx,USy;
    double Ds,Dn,De,Dw;
    double DSx,DSy;
    double denom;

         
    // reconstruct eta
    reconstruct_2D(p, d, eta, ETAs, ETAn, ETAe, ETAw);
    
    
    // reconstruct U and V
    
    LOOP
    {
    // water level       
    Ds = ETAs(i,j) + 0.5*(d->depth(i,j) + d->depth(i-1,j));
    Dn = ETAn(i,j) + 0.5*(d->depth(i,j) + d->depth(i+1,j));
    De = ETAe(i,j) + 0.5*(d->depth(i,j) + d->depth(i,j-1));
    Dw = ETAw(i,j) + 0.5*(d->depth(i,j) + d->depth(i,j+1));
    
    USx = 0.5*(Fs[IJK]+Fn[IJK]) + sqrt(9.81*Ds) - sqrt(9.81*Dn);
    USy = 0.5*(Fe[IJK]+Fw[IJK]) + sqrt(9.81*De) - sqrt(9.81*Dw);
    
    // wave speed
    Ss = MIN(Fs[IJK] - sqrt(9.81*Ds), USx - sqrt(9.81*DSx));
    Sn = MIN(Fn[IJK] + sqrt(9.81*Dn), USx + sqrt(9.81*DSx));
    
    Se = MIN(Fe[IJK] - sqrt(9.81*De), USy - sqrt(9.81*DSy));
    Sw = MIN(Fw[IJK] + sqrt(9.81*Dw), USy + sqrt(9.81*DSy));
    
    
        // final flux x-dir
        if(Ss>=0.0)
        Fx[IJK] = Fs[IJK];
        
        else
        if(Sn<=0.0)
        Fx[IJK] = Fn[IJK];
        
        else
        {
        denom = Sn-Ss;
        denom = denom>1.0e-10?denom:1.0e10;
        
        Fx[IJK] = (Sn*Fs[IJK] - Ss*Fn[IJK] + Sn*Ss*(Fn[IJK] - Fs[IJK]))/denom;
        }
        
        // final flux y-dir
        if(Se>=0.0)
        Fy[IJK] = Fe[IJK];
        
        else
        if(Sw<=0.0)
        Fy[IJK] = Fw[IJK];
        
        else
        {
        denom = Sw-Se;
        denom = denom>1.0e-10?denom:1.0e10;
        
        Fy[IJK] = (Sw*Fe[IJK] - Se*Fw[IJK] + Sw*Se*(Fw[IJK] - Fw[IJK]))/denom;
        }
    }
	
}


