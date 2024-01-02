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

#include"pressure_reference.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void pressure_reference::reference_ini(lexer*p, fdm* a, ghostcell *pgc)
{
    double xmin,xmax,ymin,zmax;
    int gcinglobal=0;
    int gcoutglobal=0;
    
    
    
    gcinglobal=pgc->globalisum(p->gcin_count);
    gcoutglobal=pgc->globalisum(p->gcout_count);
    
    //cout<<p->mpirank<<" gcinglobal: "<<gcinglobal<<endl;
    
    // ini gage location
    if(((p->B32==0 && p->B30==1) || p->B30==2 || p->B30==3))// && p->B99<3)
    {
        xmin=ymin=1.0e20;
        zmax=-1.0e20;
    
        //find active smallest xy inlet location
        if(gcinglobal>0)
        for(n=0;n<p->gcin_count;n++)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
        xmin = MIN(xmin,p->XP[IP]);
        ymin = MIN(ymin,p->YP[JP]);
        zmax = MAX(zmax,p->ZP[KP]);
        }
        
        if(gcinglobal==0)
        LOOP
        {    
        xmin = MIN(xmin,p->XP[IP]);
        ymin = MIN(ymin,p->YP[JP]);
        zmax = MAX(zmax,p->ZP[KP]);
        }
        
        p->B32_x = pgc->globalmin(xmin);
        p->B32_y = pgc->globalmin(ymin);
        p->B32_z = pgc->globalmax(zmax);
        
        if(p->mpirank==0)
        cout<<"pressure gage location  x: "<<p->B32_x<<" y: "<<p->B32_y<<" z: "<<p->B32_z<<endl;
    }
    
    // ini gage location AWA beach
    if(((p->B32==0 && p->B30==1) || p->B30==2 || p->B30==3) && p->B99>=3)
    {
        ymin=1.0e20;
        xmax=zmax=-1.0e20;
        
        //find active smallest xy inlet location
        if(gcoutglobal>0)
        for(n=0;n<p->gcout_count;n++)
        {
        i=p->gcout[n][0];
        j=p->gcout[n][1];
        k=p->gcout[n][2];
        
        xmax = MAX(xmax,p->XP[IP]);
        ymin = MIN(ymin,p->YP[JP]);
        zmax = MAX(zmax,p->ZP[KP]);
        }
        
        if(gcinglobal==0)
        LOOP
        {    
        xmax = MAX(xmax,p->XP[IP]);
        ymin = MIN(ymin,p->YP[JP]);
        zmax = MAX(zmax,p->ZP[KP]);
        }
        
        p->B32_x = pgc->globalmax(xmax);
        p->B32_y = pgc->globalmin(ymin);
        p->B32_z = pgc->globalmax(zmax);
        
        if(p->mpirank==0)
        cout<<"pressure gage location  x: "<<p->B32_x<<" y: "<<p->B32_y<<" z: "<<p->B32_z<<endl;
    }
    
    
    //atmospheric reference pressure
    if(p->B30==3)
    {
    gageval = -1.0e20;
    
    if(p->B32_x>=p->originx && p->B32_x<p->endx)
    if(p->B32_y>=p->originy && p->B32_y<p->endy)
    if(p->B32_z>=p->originz && p->B32_z<p->endz)
    gageval = p->ccipol4(a->press,p->B32_x,p->B32_y,p->B32_z);
    
    p->B31 = pgc->globalmax(gageval);
    }
    
}



