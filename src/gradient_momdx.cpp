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

#include"gradient.h"
#include"fdm.h"
#include"lexer.h"

//momentum->momentum


// *************************************************************************************
// *************************************************************************************
// X X X X X X X
// *************************************************************************************
// *************************************************************************************

// **********************************************************
// UXDX2
// **********************************************************

double gradient::udx(fdm* a)
{
	grad = (a->u(i+1,j,k) - a->u(i-1,j,k))/(p->DXN[IP]+p->DXN[IP1]);

	return grad;
}

double gradient::udy(fdm* a)
{
	grad = (a->u(i,j+1,k) - a->u(i,j-1,k))/(p->DYP[JP]+p->DXP[IM1]);

	return grad;
}

double gradient::udz(fdm* a)
{
	grad = (a->u(i,j+1,k) - a->u(i,j-1,k))/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

// **********************************************************
// UDXX2
// **********************************************************

double gradient::udxx(fdm* a)
{
	grad = ((a->u(i+1,j,k) - a->u(i,j,k))/p->DXN[IP1] - (a->u(i,j,k) + a->u(i-1,j,k))/p->DXN[IP])/p->DXP[IP];

	return grad;
}

double gradient::udyy(fdm* a)
{
    grad = ((a->u(i,j+1,k) - a->u(i,j,k))/p->DYP[JP] - (a->u(i,j,k) + a->u(i,j-1,k))/p->DYP[JM1])/p->DYN[JP];

	return grad;
}

double gradient::udzz(fdm* a)
{
	grad = ((a->u(i,j,k+1) - a->u(i,j,k))/p->DZP[KP] - (a->u(i,j,k) + a->u(i,j,k-1))/p->DZP[KM1])/p->DZN[KP];

	return grad;
}


// *************************************************************************************
// *************************************************************************************
// Y Y Y Y Y Y
// *************************************************************************************
// *************************************************************************************

// **********************************************************
// VDX2
// **********************************************************

double gradient::vdx(fdm* a)
{
	grad = (a->v(i+1,j,k) - a->v(i-1,j,k))/(p->DXP[IP]+p->DXP[IP1]);

	return grad;
}

double gradient::vdy(fdm* a)
{
	grad = (a->v(i,j+1,k) - a->v(i,j-1,k))/(p->DYN[JP]+p->DYN[JP1]);

	return grad;
}


double gradient::vdz(fdm* a)
{
	grad = (a->v(i,j+1,k) - a->v(i,j-1,k))/(p->DZP[KP]+p->DZP[KP1]);

	return grad;
}

// **********************************************************
// VDXX2
// **********************************************************

double gradient::vdxx(fdm* a)
{ 
    grad = ((a->v(i+1,j,k) - a->v(i,j,k))/p->DXP[IP] - (a->v(i,j,k) + a->v(i-1,j,k))/p->DXP[IM1])/p->DXN[IP];

	return grad;
}

double gradient::vdyy(fdm* a)
{
	grad = ((a->v(i,j+1,k) - a->v(i,j,k))/p->DYN[JP1] - (a->v(i,j,k) + a->v(i,j-1,k))/p->DYN[JP])/p->DYP[JP];

	return grad;
}


double gradient::vdzz(fdm* a)
{
	grad = ((a->v(i,j,k+1) - a->v(i,j,k))/p->DZP[KP] - (a->v(i,j,k) + a->v(i,j,k-1))/p->DZP[KM1])/p->DZN[KP];

	return grad;
}

// *************************************************************************************
// *************************************************************************************
// Z Z Z Z Z
// *************************************************************************************
// *************************************************************************************

// **********************************************************
// ZX2
// **********************************************************

double gradient::wdx(fdm* a)
{
	grad = (a->w(i+1,j,k) - a->w(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);

	return grad;
}

double gradient::wdy(fdm* a)
{
	grad = (a->w(i,j+1,k) - a->w(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]);

	return grad;
}


double gradient::wdz(fdm* a)
{
	grad = (a->w(i,j,k+1) - a->w(i,j,k-1))/(p->DZN[KP]+p->DZN[KP1]);

	return grad;
}

// **********************************************************
// WDXX2
// **********************************************************
double gradient::wdxx(fdm* a)
{
    grad = ((a->w(i+1,j,k) - a->w(i,j,k))/p->DXP[IP] - (a->w(i,j,k) + a->w(i-1,j,k))/p->DXP[IM1])/p->DXN[IP];

	return grad;
}

double gradient::wdyy(fdm* a)
{
    grad = ((a->w(i,j+1,k) - a->w(i,j,k))/p->DYP[JP] - (a->w(i,j,k) + a->w(i,j-1,k))/p->DYP[JM1])/p->DYN[JP];

	return grad;
}


double gradient::wdzz(fdm* a)
{
    grad = ((a->w(i,j,k+1) - a->w(i,j,k))/p->DZN[KP1] - (a->w(i,j,k) + a->w(i,j,k-1))/p->DZN[KP])/p->DZP[KP];

	return grad;
}


