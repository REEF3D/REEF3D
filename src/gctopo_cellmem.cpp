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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"


void ghostcell::velcell_update(lexer *p, fdm *a, int **cellmem, int cellcount,double xdelt, double ydelt, double zdelt, int dir)
{
    int nn;
	
	// solid->fluid
    for(nn=0;nn<cellcount;++nn)
    if(cellmem[nn][3]==1)
    {
    i=cellmem[nn][0];
    j=cellmem[nn][1];
    k=cellmem[nn][2];

    if(dir==1)
	a->u(i,j,k)=0.95*a->u(i,j,k+1);

	if(dir==2)
	a->v(i,j,k)=0.95*a->v(i,j,k+1);

	if(dir==3)
	a->w(i,j,k)=0.0;
    }
	
	// fluid->solid
	for(nn=0;nn<cellcount;++nn)
    if(cellmem[nn][3]==2)
    {
    i=cellmem[nn][0];
    j=cellmem[nn][1];
    k=cellmem[nn][2];

			if(dir==1)
			{			
			if(p->flag1[IJKp1]>0) 
			a->u(i,j,k+1)=0.5*a->u(i,j,k+1);
			}

			if(dir==2)
			{
			if(p->flag2[IJKp1]>0) 
			a->v(i,j,k+1)=0.5*a->v(i,j,k+1);			
			}

			if(dir==3)
			{
			if(p->flag3[IJKp1]>0) 
			a->w(i,j,k+1)=0.5*a->w(i,j,k+1);
			}
    }
    
    // NEW vertical interpolation update
    
    

}

void ghostcell::gctopo_scalarupdate(lexer *p, fdm *a, int **cellmem, int cellcount, field &f)
{
    int nn;
	double nx,ny,nz,norm;
	double posx,posy,posz;
	double locx,locy,locz;
	double topoval,fval;

	// solid->fluid
    for(nn=0;nn<cellcount;++nn)
	if(cellmem[nn][3]==1)
    {
    i=cellmem[nn][0];
    j=cellmem[nn][1];
    k=cellmem[nn][2];
	
	f(i,j,k)=0.75*f(i,j,k+1);
    }
	
	// fluid->solid
	for(nn=0;nn<cellcount;++nn)
	if(cellmem[nn][3]==2)
    {
    i=cellmem[nn][0];
    j=cellmem[nn][1];
    k=cellmem[nn][2];


	f(i,j,k+1)=0.1*f(i,j,k+1);
    }
}

void ghostcell::gctopo_pressureupdate(lexer *p, fdm *a, int **cellmem, int cellcount, field &f)
{
    int nn;
	double nx,ny,nz,norm;
	double posx,posy,posz;
	double locx,locy,locz;
	double topoval,fval;

    // solid->fluid
    for(nn=0;nn<cellcount;++nn)
	if(cellmem[nn][3]==1)
    {
    i=cellmem[nn][0];
    j=cellmem[nn][1];
    k=cellmem[nn][2];

	
	fval = a->press(i,j,k+1) + p->DZP[KP]*a->ro(i,j,k)*fabs(p->W22);

	f(i,j,k)=fval;
    }

}
