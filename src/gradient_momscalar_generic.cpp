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

#include"gradient.h"
#include"fdm.h"
#include"lexer.h"
// **********************************************************
// XDX2
// **********************************************************

double gradient::pudx(lexer *p, field &u)
{
	pip=1;
	grad = (u(i,j,k) - u(i-1,j,k))/p->DXN[IP];
	pip=0;

	return grad;
}

double gradient::pudy(lexer *p, field &u)
{
    f1 = f2 = 0.5;
    
    if(p->flag1[IJK]<0 && p->flag1[IJm1K]<0 && p->flag1[IJp1K]<0)
    {
    f2 = 0.0;
    f1 = 1.0;
    }
    
    if(p->flag1[Im1JK]<0 && p->flag1[Im1Jm1K]<0 && p->flag1[Im1Jp1K]<0)
    {
    f2 = 1.0;
    f1 = 0.0;
    }
    
	pip=2;
	grad = ((f2*u(i,j+1,k)+f1*u(i-1,j+1,k)) - (f2*u(i,j-1,k)+f1*u(i-1,j-1,k)))/(p->DYP[JP]+p->DYP[JM1]);
	pip=0;

	return grad;
}

double gradient::pudz(lexer *p, field &u)
{
    f1 = f2 = 0.5;
    
    if(p->flag1[IJK]<0 && p->flag1[IJKm1]<0 && p->flag1[IJKp1]<0)
    {
    f2 = 0.0;
    f1 = 1.0;
    }
    
    if(p->flag1[Im1JK]<0 && p->flag1[Im1JKm1]<0 && p->flag1[Im1JKp1]<0)
    {
    f2 = 1.0;
    f1 = 0.0;
    }
    
	pip=3;
	grad = ((f2*u(i,j,k+1)+f1*u(i-1,j,k+1)) - (f2*u(i,j,k-1)+f1*u(i-1,j,k-1)))/(p->DZP[KP]+p->DZP[KM1]);
	pip=0;

	return grad;
}

// **********************************************************
// YDX2
// **********************************************************

double gradient::pvdx(lexer *p, field &v)
{
    f1 = f2 = 0.5;
    
    if(p->flag2[IJK]<0 && p->flag2[Im1JK]<0 && p->flag2[Ip1JK]<0)
    {
    f2 = 0.0;
    f1 = 1.0;
    }
    
    if(p->flag2[IJm1K]<0 && p->flag2[Im1Jm1K]<0 && p->flag2[Ip1Jm1K]<0)
    {
    f2 = 1.0;
    f1 = 0.0;
    }
    
	pip=1;
	grad = ((f2*v(i+1,j,k)+f1*v(i+1,j-1,k)) - (f2*v(i-1,j,k)+f1*v(i-1,j-1,k)))/(p->DXP[IP]+p->DXP[IM1]);
	pip=0;

	return grad;
}

double gradient::pvdy(lexer *p, field &v)
{
	pip=2;
	grad = (v(i,j,k) - v(i,j-1,k))/(p->DYN[JP]);
	pip=0;

	return grad;
}

double gradient::pvdz(lexer *p, field &v)
{
    f1 = f2 = 0.5;
    
    if(p->flag2[IJK]<0 && p->flag2[IJKm1]<0 && p->flag2[IJKp1]<0)
    {
    f2 = 0.0;
    f1 = 1.0;
    }
    
    if(p->flag2[IJm1K]<0 && p->flag2[IJm1Km1]<0 && p->flag2[IJm1Kp1]<0)
    {
    f2 = 1.0;
    f1 = 0.0;
    }
    
	pip=3;
	grad = ((f2*v(i,j,k+1)+f1*v(i,j-1,k+1)) - (f2*v(i,j,k-1)+f1*v(i,j-1,k-1)))/(p->DZP[KP]+p->DZP[KM1]); 
	pip=0;

	return grad;
}

// **********************************************************
// ZX2
// **********************************************************

double gradient::pwdx(lexer *p, field &w)
{
    f1 = f2 = 0.5;
    
    if(p->flag3[IJK]<0 && p->flag3[Im1JK]<0 && p->flag3[Ip1JK]<0)
    {
    f2 = 0.0;
    f1 = 1.0;
    }
    
    if(p->flag3[IJKm1]<0 && p->flag3[Im1JKm1]<0 && p->flag3[Ip1JKm1]<0)
    {
    f2 = 1.0;
    f1 = 0.0;
    }
    
	pip=1;
	grad = ((f2*w(i+1,j,k)+f1*w(i+1,j,k-1)) - (f2*w(i-1,j,k)+f1*w(i-1,j,k-1)))/(p->DXP[IP]+p->DXP[IM1]);
	pip=0;

	return grad;
}

double gradient::pwdy(lexer *p, field &w)
{
    f1 = f2 = 0.5;
    
    if(p->flag3[IJK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJp1K]<0)
    {
    f2 = 0.0;
    f1 = 1.0;
    }
    
    if(p->flag3[IJKm1]<0 && p->flag3[IJm1Km1]<0 && p->flag3[IJp1Km1]<0)
    {
    f2 = 1.0;
    f1 = 0.0;
    }
    
	pip=2;
	grad = ((f2*w(i,j+1,k)+f1*w(i,j+1,k-1)) - (f2*w(i,j-1,k)+f1*w(i,j-1,k-1)))/(p->DYP[JP]+p->DYP[JM1]);
	pip=0;

	return grad;
}

double gradient::pwdz(lexer *p, field &w)
{
	pip=3;
	grad = (w(i,j,k) - w(i,j,k-1))/(p->DZN[KP]);
	pip=0;

	return grad;
}

