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

#include"patchBC.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

void patchBC::patchBC_pressure(lexer *p, fdm *a, ghostcell *pgc, field &press)
{

        
    for(qq=0;qq<obj_count;++qq)
    if(patch[qq]->pressure_flag==1)
    for(n=0;n<patch[qq]->gcb_count;++n)
    {
    i=patch[qq]->gcb[n][0];
    j=patch[qq]->gcb[n][1];
    k=patch[qq]->gcb[n][2];
    
    double eps,H;
                
    eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
    if(a->phi(i,j,k)>eps)
    H=1.0;

    if(a->phi(i,j,k)<-eps)
    H=0.0;

    if(fabs(a->phi(i,j,k))<=eps)
    H=0.5*(1.0 + a->phi(i,j,k)/eps + (1.0/PI)*sin((PI*a->phi(i,j,k))/eps));
        
    //pval=(1.0-H)*a->press(i,j,k);
    
        if(patch[qq]->gcb[n][3]==1)
        {
        press(i,j,k)   =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i-1,j,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i-2,j,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i-3,j,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        press(i,j,k)   =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j+1,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j+2,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j+3,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        press(i,j,k)   =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j-1,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j-2,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j-3,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        press(i,j,k)   =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i+1,j,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i+2,j,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i+3,j,k) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==5)
        {
        press(i,j,k)   =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j,k-1) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j,k-2) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j,k-3) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        }
        
        if(patch[qq]->gcb[n][3]==6)
        {
        press(i,j,k)   =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j,k+1) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j,k+2) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        press(i,j,k+3) =  H*patch[qq]->pressure + (1.0-H)*press(i,j,k);
        }
    
    }
} 

void patchBC::patchBC_pressure2D(lexer*, ghostcell*, slice&)
{
}

void patchBC::patchBC_pressure2D_ugrad(lexer*, fdm2D*, slice&,slice&)
{
}

void patchBC::patchBC_pressure2D_vgrad(lexer*, fdm2D*, slice&, slice&)
{
}