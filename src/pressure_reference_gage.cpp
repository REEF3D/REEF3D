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

#include"pressure_reference.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void pressure_reference::gage_fixed(lexer*p, fdm* a, ghostcell *pgc)
{
    gageval = -1.0e20;
    
    if(p->B32_x>=p->originx && p->B32_x<p->endx)
    if(p->B32_y>=p->originy && p->B32_y<p->endy)
    if(p->B32_z>=p->originz && p->B32_z<p->endz)
    gageval = p->ccipol4(a->press,p->B32_x,p->B32_y,p->B32_z);
    
    gageval = pgc->globalmax(gageval);
    
    if(p->B33==1)
    p->pressgage = gageval - p->B31;
    
    if(p->B33==2 && gageval>-0.9e20)
    LOOP
    a->press(i,j,k) -= (gageval - p->B31);
    
    //cout<<p->mpirank<<" gage: "<<gageval<<" x: "<<p->B32_x<<" y: "<<p->B32_y<<" z: "<<p->B32_z<<endl;
    
}

void pressure_reference::gage_fsf(lexer*p, fdm* a, ghostcell *pgc)
{
    p->B32_z = gageval = -1.0e20;
    
    if(p->B32_x>=p->originx && p->B32_x<p->endx)
    if(p->B32_y>=p->originy && p->B32_y<p->endy)
    {
    i =  p->posc_i(p->B32_x);
    j =  p->posc_j(p->B32_y);
    
        for(k=0; k<p->knoz-1; ++k)
        PCHECK
        if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
        p->B32_z = MAX(p->B32_z,-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k) + 1.0e-8) + p->pos_z());
        
        
    //cout<<p->mpirank<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" z: "<<p->B32_z<<endl;
        
    gageval = p->ccipol4(a->press,p->B32_x,p->B32_y,p->B32_z);
    }
    
    p->B32_z = pgc->globalmax(p->B32_z);
    
    gageval = pgc->globalmax(gageval);
    
    if(p->B33==1)
    p->pressgage = gageval - p->B31;
    
    if(p->B33==2 && gageval>-0.9e20)
    LOOP
    a->press(i,j,k) -= (gageval - p->B31);

    if(p->mpirank==0)
    cout<<" gage_fsf: "<<gageval<<" x: "<<p->B32_x<<" y: "<<p->B32_y<<" z: "<<p->B32_z<<endl;
    
}