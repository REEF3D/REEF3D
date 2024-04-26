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
#include"lexer.h"
#include"field.h"
#include"vec.h"
#include"fdm.h"

void ghostcell::neumann(field& f, int gcv, int bc, int cs)
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

void ghostcell::gcV_neumann(vec &x, int gcv, int bc, int cs, int id)
{
    n=id;
    
	if(cs==1)
    {
	x.V[Im1_J_K_4]=x.V[I_J_K_4];
    x.V[Im2_J_K_4]=x.V[I_J_K_4];
    x.V[Im3_J_K_4]=x.V[I_J_K_4];
    }

	if(cs==2)
    {
	x.V[I_Jp1_K_4]=x.V[I_J_K_4];
    x.V[I_Jp2_K_4]=x.V[I_J_K_4];
    x.V[I_Jp3_K_4]=x.V[I_J_K_4];
    }
    
	if(cs==3)
	{
	x.V[I_Jm1_K_4]=x.V[I_J_K_4];
    x.V[I_Jm2_K_4]=x.V[I_J_K_4];
    x.V[I_Jm3_K_4]=x.V[I_J_K_4];
    }

	if(cs==4)
	{
    x.V[Ip1_J_K_4]=x.V[I_J_K_4];
    x.V[Ip2_J_K_4]=x.V[I_J_K_4];
    x.V[Ip3_J_K_4]=x.V[I_J_K_4];
    }

	if(cs==5)
	{
    x.V[I_J_Km1_4]=x.V[I_J_K_4];
    x.V[I_J_Km2_4]=x.V[I_J_K_4];
    x.V[I_J_Km3_4]=x.V[I_J_K_4];
    }

	if(cs==6)
	{
    x.V[I_J_Kp1_4]=x.V[I_J_K_4];
    x.V[I_J_Kp2_4]=x.V[I_J_K_4];
    x.V[I_J_Kp3_4]=x.V[I_J_K_4];
    }
}

void ghostcell::neumannV(double *f, int gcv, int bc, int cs)
{
	if(cs==1)
    {
	f[Im1JK]=f[IJK];
    f[Im2JK]=f[IJK];
    f[Im2JK]=f[IJK];
    }
    
    if(cs==2)
    {
	f[IJp1K]=f[IJK];
    f[IJp2K]=f[IJK];
    f[IJp3K]=f[IJK];
    }
    
    if(cs==3)
    {
	f[IJm1K]=f[IJK];
    f[IJm2K]=f[IJK];
    f[IJm3K]=f[IJK];
    
    cout<<"NEUMAN"<<endl;
    }
    
    if(cs==4)
    {
	f[Ip1JK]=f[IJK];
    f[Ip2JK]=f[IJK];
    f[Ip3JK]=f[IJK];
    }
    
    if(cs==5)
    {
	f[IJKm1]=f[IJK];
    f[IJKm2]=f[IJK];
    f[IJKm3]=f[IJK];
    }
    
    if(cs==6)
    {
	f[IJKp1]=f[IJK];
    f[IJKp2]=f[IJK];
    f[IJKp3]=f[IJK];
    }
}
