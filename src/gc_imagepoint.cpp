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
#include"field.h"
#include"lexer.h"

void ghostcell::imagepoint(lexer *p, field& f,double &x_ip, double& val_ip, double dist, int cs)
{
    double dx;
    
    if(cs==1||cs==4)
    dx = p->DXP[IP];
    
    if(cs==2||cs==3)
    dx = p->DYP[JP];
    
    if(cs==5||cs==6)
    dx = p->DZP[KP];
    
	double y0,y1;
	y1=0.0;      // x_j-1
    y0=f(i,j,k); // x_j

	//fill y[]
	if(cs==1 )
	y1=f(i+1,j,k);

	if(cs==2)
	y1=f(i,j-1,k);

	if(cs==3)
	y1=f(i,j+1,k);

	if(cs==4)
	y1=f(i-1,j,k);

	if(cs==5)
	y1=f(i,j,k+1);

	if(cs==6)
	y1=f(i,j,k-1);


	x_ip = -(gamma*dx);

	val_ip= (1.0-gamma)*y0  + gamma*y1;
}

