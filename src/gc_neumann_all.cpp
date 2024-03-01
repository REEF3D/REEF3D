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

#include"ghostcell.h"
#include"field.h"
#include"fdm.h"
#include<math.h>

void ghostcell::neumann_all(field& f, int gcv, int bc, int cs)
{
	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=f(i,j,k);

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=f(i,j,k);

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=f(i,j,k);

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=f(i,j,k);

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=f(i,j,k);

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=f(i,j,k);
}

void ghostcell::gcV_neumann_all(vec &x, int gcv, int bc, int cs, int id)
{
	n=id;
    
	if(cs==1)
    {
	x.V[Im1_J_K_4a]=x.V[I_J_K_4a];
    x.V[Im2_J_K_4a]=x.V[I_J_K_4a];
    x.V[Im3_J_K_4a]=x.V[I_J_K_4a];
    }

	if(cs==2)
    {
	x.V[I_Jp1_K_4a]=x.V[I_J_K_4a];
    x.V[I_Jp2_K_4a]=x.V[I_J_K_4a];
    x.V[I_Jp3_K_4a]=x.V[I_J_K_4a];
    }
    
	if(cs==3)
	{
	x.V[I_Jm1_K_4a]=x.V[I_J_K_4a];
    x.V[I_Jm2_K_4a]=x.V[I_J_K_4a];
    x.V[I_Jm3_K_4a]=x.V[I_J_K_4a];
    }

	if(cs==4)
	{
    x.V[Ip1_J_K_4a]=x.V[I_J_K_4a];
    x.V[Ip2_J_K_4a]=x.V[I_J_K_4a];
    x.V[Ip3_J_K_4a]=x.V[I_J_K_4a];
    }

	if(cs==5)
	{
    x.V[I_J_Km1_4a]=x.V[I_J_K_4a];
    x.V[I_J_Km2_4a]=x.V[I_J_K_4a];
    x.V[I_J_Km3_4a]=x.V[I_J_K_4a];
    }

	if(cs==6)
	{
    x.V[I_J_Kp1_4a]=x.V[I_J_K_4a];
    x.V[I_J_Kp2_4a]=x.V[I_J_K_4a];
    x.V[I_J_Kp3_4a]=x.V[I_J_K_4a];
    }
}

void ghostcell::gcV_neumann_6V(vec &x, int gcv, int bc, int cs, int id)
{
	n=id;
    
	if(cs==1)
    {
	x.V[Im1_J_K_6]=x.V[I_J_K_6];
    x.V[Im2_J_K_6]=x.V[I_J_K_6];
    x.V[Im3_J_K_6]=x.V[I_J_K_6];
    }

	if(cs==2)
    {
	x.V[I_Jp1_K_6]=x.V[I_J_K_6];
    x.V[I_Jp2_K_6]=x.V[I_J_K_6];
    x.V[I_Jp3_K_6]=x.V[I_J_K_6];
    }
    
	if(cs==3)
	{
	x.V[I_Jm1_K_6]=x.V[I_J_K_6];
    x.V[I_Jm2_K_6]=x.V[I_J_K_6];
    x.V[I_Jm3_K_6]=x.V[I_J_K_6];
    }

	if(cs==4)
	{
    x.V[Ip1_J_K_6]=x.V[I_J_K_6];
    x.V[Ip2_J_K_6]=x.V[I_J_K_6];
    x.V[Ip3_J_K_6]=x.V[I_J_K_6];
    }

	if(cs==5)
	{
    x.V[I_J_Km1_6]=x.V[I_J_K_6];
    x.V[I_J_Km2_6]=x.V[I_J_K_6];
    x.V[I_J_Km3_6]=x.V[I_J_K_6];
    }

	if(cs==6)
	{
    x.V[I_J_Kp1_6]=x.V[I_J_K_6];
    x.V[I_J_Kp2_6]=x.V[I_J_K_6];
    x.V[I_J_Kp3_6]=x.V[I_J_K_6];
    }
}
