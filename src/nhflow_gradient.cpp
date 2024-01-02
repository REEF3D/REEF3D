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


#include"nhflow_gradient.h"
#include"fdm_nhf.h"
#include"lexer.h"

nhflow_gradient::nhflow_gradient(lexer* pp) : tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(0.000001),smallnum(1.0e-20),dx(pp->DXM)
{
    p=pp;
    
    grad=0.0;
}

nhflow_gradient::~nhflow_gradient()
{

}

// **********************************************************
// DUXDX2
// **********************************************************

double nhflow_gradient::dudx(double *U)
{
	//grad = (a->u(i+1,j,k) - a->u(i-1,j,k))/(p->DXN[IP]+p->DXN[IP1]);

	return grad;
}

double nhflow_gradient::dudy(double *U)
{
	//grad = (a->u(i,j+1,k) - a->u(i,j-1,k))/(p->DYP[JP]+p->DXP[IM1]);

	return grad;
}

double nhflow_gradient::dudz(double *U)
{
	//grad = (a->u(i,j+1,k) - a->u(i,j-1,k))/(p->DZP[KP]+p->DZP[KM1]);

	return grad;
}

// **********************************************************
// DUDXX2
// **********************************************************

double nhflow_gradient::dudxx(double *U)
{
	//grad = ((a->u(i+1,j,k) - a->u(i,j,k))/p->DXN[IP1] - (a->u(i,j,k) + a->u(i-1,j,k))/p->DXN[IP])/p->DXP[IP];

	return grad;
}

double nhflow_gradient::dudyy(double *U)
{
    //grad = ((a->u(i,j+1,k) - a->u(i,j,k))/p->DYP[JP] - (a->u(i,j,k) + a->u(i,j-1,k))/p->DYP[JM1])/p->DYN[JP];

	return grad;
}

double nhflow_gradient::dudzz(double *U)
{
	//grad = ((a->u(i,j,k+1) - a->u(i,j,k))/p->DZP[KP] - (a->u(i,j,k) + a->u(i,j,k-1))/p->DZP[KM1])/p->DZN[KP];

	return grad;
}

// **********************************************************
// DVDX2
// **********************************************************

double nhflow_gradient::dvdx(double *V)
{
	//grad = (a->v(i+1,j,k) - a->v(i-1,j,k))/(p->DXP[IP]+p->DXP[IP1]);

	return grad;
}

double nhflow_gradient::dvdy(double *V)
{
	//grad = (a->v(i,j+1,k) - a->v(i,j-1,k))/(p->DYN[JP]+p->DYN[JP1]);

	return grad;
}


double nhflow_gradient::dvdz(double *V)
{
	//grad = (a->v(i,j+1,k) - a->v(i,j-1,k))/(p->DZP[KP]+p->DZP[KP1]);

	return grad;
}

// **********************************************************
// DVDXX2
// **********************************************************

double nhflow_gradient::dvdxx(double *V)
{ 
    //grad = ((a->v(i+1,j,k) - a->v(i,j,k))/p->DXP[IP] - (a->v(i,j,k) + a->v(i-1,j,k))/p->DXP[IM1])/p->DXN[IP];

	return grad;
}

double nhflow_gradient::dvdyy(double *V)
{
	//grad = ((a->v(i,j+1,k) - a->v(i,j,k))/p->DYN[JP1] - (a->v(i,j,k) + a->v(i,j-1,k))/p->DYN[JP])/p->DYP[JP];

	return grad;
}


double nhflow_gradient::dvdzz(double *V)
{
	//grad = ((a->v(i,j,k+1) - a->v(i,j,k))/p->DZP[KP] - (a->v(i,j,k) + a->v(i,j,k-1))/p->DZP[KM1])/p->DZN[KP];

	return grad;
}

// **********************************************************
// DZX2
// **********************************************************

double nhflow_gradient::dwdx(double *W)
{
	//grad = (a->w(i+1,j,k) - a->w(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);

	return grad;
}

double nhflow_gradient::dwdy(double *W)
{
	//grad = (a->w(i,j+1,k) - a->w(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]);

	return grad;
}


double nhflow_gradient::dwdz(double *W)
{
	//grad = (a->w(i,j,k+1) - a->w(i,j,k-1))/(p->DZN[KP]+p->DZN[KP1]);

	return grad;
}

// **********************************************************
// DWDXX2
// **********************************************************
double nhflow_gradient::dwdxx(double *W)
{
    //grad = ((a->w(i+1,j,k) - a->w(i,j,k))/p->DXP[IP] - (a->w(i,j,k) + a->w(i-1,j,k))/p->DXP[IM1])/p->DXN[IP];

	return grad;
}

double nhflow_gradient::dwdyy(double *W)
{
    //grad = ((a->w(i,j+1,k) - a->w(i,j,k))/p->DYP[JP] - (a->w(i,j,k) + a->w(i,j-1,k))/p->DYP[JM1])/p->DYN[JP];

	return grad;
}


double nhflow_gradient::dwdzz(double *W)
{
    //grad = ((a->w(i,j,k+1) - a->w(i,j,k))/p->DZN[KP1] - (a->w(i,j,k) + a->w(i,j,k-1))/p->DZN[KP])/p->DZP[KP];

	return grad;
}


