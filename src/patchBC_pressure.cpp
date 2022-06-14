/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
    
        if(patch[qq]->gcb[n][3]==1)
        {
        press(i-1,j,k) =  patch[qq]->pressure;
        press(i-2,j,k) =  patch[qq]->pressure;
        press(i-3,j,k) =  patch[qq]->pressure;
        }
        
        if(patch[qq]->gcb[n][3]==2)
        {
        press(i,j+1,k) =  patch[qq]->pressure;
        press(i,j+2,k) =  patch[qq]->pressure;
        press(i,j+3,k) =  patch[qq]->pressure;
        }
        
        if(patch[qq]->gcb[n][3]==3)
        {
        press(i,j-1,k) =  patch[qq]->pressure;
        press(i,j-2,k) =  patch[qq]->pressure;
        press(i,j-3,k) =  patch[qq]->pressure;
        }
        
        if(patch[qq]->gcb[n][3]==4)
        {
        press(i+1,j,k) =  patch[qq]->pressure;
        press(i+2,j,k) =  patch[qq]->pressure;
        press(i+3,j,k) =  patch[qq]->pressure;
        }
        
        if(patch[qq]->gcb[n][3]==5)
        {
        //cout<<p->mpirank<<" TEST PRESS"<<endl;
        press(i,j,k-1) =  patch[qq]->pressure;
        press(i,j,k-2) =  patch[qq]->pressure;
        press(i,j,k-3) =  patch[qq]->pressure;
        }
        
        if(patch[qq]->gcb[n][3]==6)
        {
        press(i,j,k+1) =  patch[qq]->pressure;
        press(i,j,k+2) =  patch[qq]->pressure;
        press(i,j,k+3) =  patch[qq]->pressure;
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