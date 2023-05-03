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

#include"lexer.h"
#include"ghostcell.h"
#include"field.h"
#include"fdm.h"

void ghostcell::lsm(lexer *p,field& f,double dist,int gcv, int bc, int cs)
{
	if(cs==1)
	for(q=1;q<=margin;++q)
	f(i-q,j,k) = double(q+1)*f(i,j,k) - double(q)*f(i+1,j,k);

	if(cs==2)
	for(q=1;q<=margin;++q)
	f(i,j+q,k) = double(q+1)*f(i,j,k) - double(q)*f(i,j-1,k);

	if(cs==3)
	for(q=1;q<=margin;++q)
	f(i,j-q,k) = double(q+1)*f(i,j,k) - double(q)*f(i,j+1,k);

	if(cs==4)
	for(q=1;q<=margin;++q)
	f(i+q,j,k) = double(q+1)*f(i,j,k) - double(q)*f(i-1,j,k);

	if(cs==5)
	for(q=1;q<=margin;++q)
	f(i,j,k-q) = double(q+1)*f(i,j,k) - double(q)*f(i,j,k+1);

	if(cs==6)
	for(q=1;q<=margin;++q)
	f(i,j,k+q) = double(q+1)*f(i,j,k) - double(q)*f(i,j,k-1);
}


void ghostcell::gcV_lsm(lexer *p,vec &x, double dist,int gcv, int bc, int cs, int id)
{
    n=id;
    
	if(cs==1)
    {
    x.V[Im1_J_K_4] = 2.0*x.V[I_J_K_4] - 1.0*x.V[Ip1_J_K_4];
    x.V[Im2_J_K_4] = 3.0*x.V[I_J_K_4] - 2.0*x.V[Ip1_J_K_4];
    x.V[Im3_J_K_4] = 4.0*x.V[I_J_K_4] - 3.0*x.V[Ip1_J_K_4];
    }
    
	if(cs==2)
	{
    x.V[I_Jp1_K_4] = 2.0*x.V[I_J_K_4] - 1.0*x.V[I_Jm1_K_4];
    x.V[I_Jp2_K_4] = 3.0*x.V[I_J_K_4] - 2.0*x.V[I_Jm1_K_4];
    x.V[I_Jp3_K_4] = 4.0*x.V[I_J_K_4] - 3.0*x.V[I_Jm1_K_4];
    }

	if(cs==3)
	{
    x.V[I_Jm1_K_4] = 2.0*x.V[I_J_K_4] - 1.0*x.V[I_Jp1_K_4];
    x.V[I_Jm2_K_4] = 3.0*x.V[I_J_K_4] - 2.0*x.V[I_Jp1_K_4];
    x.V[I_Jm3_K_4] = 4.0*x.V[I_J_K_4] - 3.0*x.V[I_Jp1_K_4];
    }

	if(cs==4)
	{
    x.V[Ip1_J_K_4] = 2.0*x.V[I_J_K_4] - 1.0*x.V[Im1_J_K_4];
    x.V[Ip2_J_K_4] = 3.0*x.V[I_J_K_4] - 2.0*x.V[Im1_J_K_4];
    x.V[Ip3_J_K_4] = 4.0*x.V[I_J_K_4] - 3.0*x.V[Im1_J_K_4];
    }

	if(cs==5)
	{
    x.V[I_J_Km1_4] = 2.0*x.V[I_J_K_4] - 1.0*x.V[I_J_Kp1_4];
    x.V[I_J_Km2_4] = 3.0*x.V[I_J_K_4] - 2.0*x.V[I_J_Kp1_4];
    x.V[I_J_Km3_4] = 4.0*x.V[I_J_K_4] - 3.0*x.V[I_J_Kp1_4];
    }

	if(cs==6)
	{
    x.V[I_J_Kp1_4] = 2.0*x.V[I_J_K_4] - 1.0*x.V[I_J_Km1_4];
    x.V[I_J_Kp2_4] = 3.0*x.V[I_J_K_4] - 2.0*x.V[I_J_Km1_4];
    x.V[I_J_Kp3_4] = 4.0*x.V[I_J_K_4] - 3.0*x.V[I_J_Km1_4];
    }
	
}
