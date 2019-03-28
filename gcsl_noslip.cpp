/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fdm2D.h"
#include"ghostcell.h"
#include"field.h"
#include"slice.h"
#include"vec2D.h"

void ghostcell::gcsl_noslip(slice &f, int gcv, int bc, int cs)
{
	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j)=0.0;

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1)=0.0;

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1)=0.0;

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j)=0.0;
}

void ghostcell::gcsl_noslipV(vec2D &x, int gcv, int bc, int cs, int id)
{
	n=id;
    
	if(cs==1)
    {
	x.V[Im1_J_4]=0.0;
    x.V[Im2_J_4]=0.0;
    x.V[Im3_J_4]=0.0;
    }

	if(cs==2)
	{
	x.V[I_Jp1_4]=0.0;
    x.V[I_Jp2_4]=0.0;
    x.V[I_Jp3_4]=0.0;
    }

	if(cs==3)
	{
	x.V[I_Jm1_4]=0.0;
    x.V[I_Jm2_4]=0.0;
    x.V[I_Jm3_4]=0.0;
    }

	if(cs==4)
	{
	x.V[Ip1_J_4]=0.0;
    x.V[Ip2_J_4]=0.0;
    x.V[Ip3_J_4]=0.0;
    }
}
