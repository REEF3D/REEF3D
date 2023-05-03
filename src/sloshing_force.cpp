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

#include"sloshing_force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>

sloshing_force::sloshing_force(lexer *p, fdm* a, ghostcell *pgc)
{
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_Force",0777);
	
    if(p->mpirank==0 && p->P101>0)
    {
    // open file
	if(p->P14==0)
    result.open("REEF3D-CFD-Sloshing-Force.dat");
	
	if(p->P14==1)
	result.open("./REEF3D_CFD_Force/REEF3D-CFD-Sloshing-Force.dat");

    result<<"time \t Fx_l \t Fx_r \t Fz \t M "<<endl;
    }

}

sloshing_force::~sloshing_force()
{
    result.close();
}

void sloshing_force::start(lexer *p, fdm *a, ghostcell *pgc)
{
    
    force(p,a,pgc);
    
    // write to file
    if(p->mpirank==0)
    result<<p->simtime<<"\t "<<Fx_l<<"\t "<<Fx_r<<"\t "<<Fz<<"\t "<<M<<endl;
}

void sloshing_force::force(lexer *p, fdm *a, ghostcell *pgc)
{
    double dist_x,dist_z;
    
    Fx_l=Fx_r=Fz=M=0.0;
    
    GC4LOOP
    if(p->gcb4[n][4]==21 || p->gcb4[n][4]==22)
    {
        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];
        
        dist_x = p->B192_3 - p->pos_x();
        dist_z = p->pos_z() - p->B192_4;
        
        if(p->gcb4[n][3]==1)
        {
        if(i+p->origin_i==0)
        Fx_l-=p->DXM*p->DXM*a->press(i,j,k);
        M+=p->DXM*p->DXM*a->press(i,j,k)*dist_z;
        }
        
        if(p->gcb4[n][3]==4)
        {
        if(i+p->origin_i==p->gknox-1)
        Fx_r+=p->DXM*p->DXM*a->press(i,j,k);
        M-=p->DXM*p->DXM*a->press(i,j,k)*dist_z;
        }
        
        if(p->gcb4[n][3]==5)
        {
        if(+p->origin_k==0)
        Fz+=p->DXM*p->DXM*a->press(i,j,k);
        M+=p->DXM*p->DXM*a->press(i,j,k)*dist_x;
        }
        
    }
    
    
    Fx_l = pgc->globalsum(Fx_l);
    Fx_r = pgc->globalsum(Fx_r);
    Fz = pgc->globalsum(Fz);
    M = pgc->globalsum(M);
}

