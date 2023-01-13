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

#include"ptf_fsfbc.h"
#include"lexer.h"
#include"fdm.h"

double ptf_fsfbc::fz(lexer *p, fdm *a, field &f, slice &Fifsf)
{
    grad=0.0;
    teta=0.0;
    
    teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
    
    //cout<<"TETA: "<<teta<<" p->ZP[KP]: "<<p->ZP[KP]<<" p->ZP[KP]+teta*p->DZN[KP]: "<<p->ZP[KP]+teta*p->DZN[KP]<<endl;
    
    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3] && p->flag4[IJKm4]>0 && p->flag4[IJKm5])
    {
        if(i+p->origin_i>0)
        grad = (-(49.0/20.0)*Fifsf(i,j) + 6.0*f(i,j,k) - 7.5*f(i,j,k-1) + (20.0/3.0)*f(i,j,k-2) - (15.0/4.0)*f(i,j,k-3) + (6.0/5.0)*f(i,j,k-4) - (1.0/6.0)*f(i,j,k-5))
          /(-(49.0/20.0)*(p->ZP[KP]+teta*p->DZN[KP]) + 6.0*p->ZP[KP] - 7.5*p->ZP[KM1] + (20.0/3.0)*p->ZP[KM2] - (15.0/4.0)*p->ZP[KM3] + (6.0/5.0)*p->ZP[KM4] - (1.0/6.0)*p->ZP[KM5]);
              
        if(i+p->origin_i==0)
        grad = (-(49.0/20.0)*f(i,j,k) + 6.0*f(i,j,k-1) - 7.5*f(i,j,k-2) + (20.0/3.0)*f(i,j,k-3) - (15.0/4.0)*f(i,j,k-4) + (6.0/5.0)*f(i,j,k-5) - (1.0/6.0)*f(i,j,k-6))
          /(-(49.0/20.0)*p->ZP[KP] + 6.0*p->ZP[KM1] - 7.5*p->ZP[KM2] + (20.0/3.0)*p->ZP[KM3] - (15.0/4.0)*p->ZP[KM4] + (6.0/5.0)*p->ZP[KM5] - (1.0/6.0)*p->ZP[KM6]);
              
        return grad;
    }
    
    //else
    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0 && p->flag4[IJKm2]>0 && p->flag4[IJKm3]>0)
    {
        if(i+p->origin_i>0)
        grad = (-(25.0/12.0)*Fifsf(i,j) + 4.0*f(i,j,k) - 3.0*f(i,j,k-1) + (4.0/3.0)*f(i,j,k-2) - 0.25*f(i,j,k-3))
              /(-(25.0/12.0)*(p->ZP[KP]+teta*p->DZN[KP]) + 4.0*p->ZP[KP] - 3.0*p->ZP[KM1] + (4.0/3.0)*p->ZP[KM2] - 0.25*p->ZP[KM3]);
              
        if(i+p->origin_i==0)
        grad = (-(25.0/12.0)*f(i,j,k) + 4.0*f(i,j,k-1) - 3.0*f(i,j,k-2) + (4.0/3.0)*f(i,j,k-3) - 0.25*f(i,j,k-4))
              /(-(25.0/12.0)*p->ZP[KP] + 4.0*p->ZP[KM1] - 3.0*p->ZP[KM2] + (4.0/3.0)*p->ZP[KM3] - 0.25*p->ZP[KM4]);
                  
        return grad;
    }
    
    else
    if(p->flag4[IJK]>0 && p->flag4[IJKm1]>0)
    {
        if(i+p->origin_i>0)
        grad = (-1.5*Fifsf(i,j) + 2.0*f(i,j,k) - 0.5*f(i,j,k-1))/(-1.5*(p->ZP[KP]+teta*p->DZN[KP]) + 2.0*p->ZP[KP] - 0.5*p->ZP[KM1]);
              
        if(i+p->origin_i==0)
        grad = (-1.5*f(i,j,k) + 2.0*f(i,j,k-1) - 0.5*f(i,j,k-2))/(-1.5*p->ZP[KP] + 2.0*p->ZP[KM1] - 0.5*p->ZP[KM2]);
             
        return grad;
    }
    
    else
    {
        if(i+p->origin_i>0)
        grad = (f(i,j,k+1) - f(i,j,k))/(p->ZP[KP]+p->DZN[KP]);
              
        if(i+p->origin_i==0)
        grad = (f(i,j,k) - f(i,j,k-1))/(p->ZP[KM1]);
            
        return grad;
    }
}


