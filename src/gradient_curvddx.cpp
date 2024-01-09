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

//scalar->momentum

// **********************************************************
// XDX2
// **********************************************************

double gradient::xdx(fdm* a, field& f)
{
    pip=1;
	grad = (f(i+1,j,k) - f(i,j,k))/(dx);
    pip=0;
	return grad;
}

double gradient::xdy(fdm* a, field& f)
{
    pip=2;
	grad = (0.5*(f(i,j+1,k)+f(i+1,j+1,k)) - 0.5*(f(i,j-1,k)+f(i+1,j-1,k)))/(2.0*dx);
    pip=0;

	return grad;
}

double gradient::xdz(fdm* a, field& f)
{
    pip=3;
	grad = (0.5*(f(i,j,k+1)+f(i+1,j,k+1)) - 0.5*(f(i,j,k-1)+f(i+1,j,k-1)))/(2.0*dx);
    pip=0;

	return grad;
}

// **********************************************************
// XDXX2
// **********************************************************

double gradient::xdxx(fdm* a, field& f)
{
    pip=1;
	grad = (f(i+2,j,k)-f(i+1,j,k)-f(i,j,k)+f(i-1,j,k))/(2.0*dx*dx);
    pip=0;

	return grad;
}

double gradient::xdyy(fdm* a, field& f)
{
    pip=2;
	grad = (0.5*(f(i,j-1,k)+f(i+1,j-1,k)) - 2.0*0.5*(f(i,j,k)+f(i+1,j,k)) + 0.5*(f(i,j+1,k)+f(i+1,j+1,k)))/(dx*dx);
    pip=0;

	return grad;
}

double gradient::xdzz(fdm* a, field& f)
{
    pip=3;
	grad = (0.5*(f(i,j,k-1)+f(i+1,j,k-1)) - 2.0*0.5*(f(i,j,k)+f(i+1,j,k)) + 0.5*(f(i,j,k+1)+f(i+1,j,k+1)))/(dx*dx);
    pip=0;

	return grad;
}

// **********************************************************
// XDXY2
// **********************************************************

double gradient::xdxy(fdm* a, field& f)
{
    pip=1;
	grad = (f(i+1,j+1,k)-f(i+1,j-1,k)-f(i,j+1,k)+f(i,j-1,k))/(2.0*dx*dx);
    pip=0;

	return grad;
}

double gradient::xdxz(fdm* a, field& f)
{
    pip=1;
	grad = (f(i+1,j,k+1)-f(i+1,j,k-1)-f(i,j,k+1)+f(i,j,k-1))/(2.0*dx*dx);
    pip=0;

	return grad;
}

double gradient::xdyz(fdm* a, field& f)
{
    pip=2;
    grad = (0.5*(f(i,j+1,k+1)+f(i+1,j+1,k+1)) - 0.5*(f(i,j-1,k+1)+f(i+1,j-1,k+1))

         -0.5*(f(i,j+1,k-1)+f(i+1,j+1,k-1)) + 0.5*(f(i,j-1,k-1)+f(i+1,j-1,k-1)))/(4.0*dx*dx);
    pip=0;

	return grad;
}

// **********************************************************
// YDX2
// **********************************************************

double gradient::ydx(fdm* a, field& f)
{
    pip=1;
	grad = (0.5*(f(i+1,j,k)+f(i+1,j+1,k)) - 0.5*(f(i-1,j,k)+f(i-1,j+1,k)))/(2.0*dx);
    pip=0;

	return grad;
}

double gradient::ydy(fdm* a, field& f)
{
    pip=2;
	grad = (f(i,j+1,k) - f(i,j,k))/(dx);
    pip=0;

	return grad;
}

double gradient::ydz(fdm* a, field& f)
{
    pip=3;
	grad = (0.5*(f(i,j,k+1)+f(i,j+1,k+1)) - 0.5*(f(i,j,k-1)+f(i,j+1,k-1)))/(2.0*dx);
    pip=0;

	return grad;
}

// **********************************************************
// YDXX2
// **********************************************************

double gradient::ydxx(fdm* a, field& f)
{
    pip=1;
	grad = (0.5*(f(i-1,j,k)+f(i-1,j+1,k)) - 2.0*0.5*(f(i,j,k)+f(i,j+1,k)) + 0.5*(f(i+1,j,k)+f(i+1,j+1,k)))/(dx*dx);
    pip=0;

	return grad;
}

double gradient::ydyy(fdm* a, field& f)
{
    pip=2;
	grad = (f(i,j+2,k)-f(i,j+1,k)-f(i,j,k)+f(i,j-1,k))/(2.0*dx*dx);
    pip=0;

	return grad;
}

double gradient::ydzz(fdm* a, field& f)
{
    pip=3;
	grad = (0.5*(f(i,j,k-1)+f(i,j+1,k-1)) - 2.0*0.5*(f(i,j,k)+f(i,j+1,k)) + 0.5*(f(i,j,k+1)+f(i,j+1,k+1)))/(dx*dx);
    pip=0;

	return grad;
}

// **********************************************************
// YDXY2
// **********************************************************

double gradient::ydxy(fdm* a, field& f)
{
    pip=2;
	grad = (f(i+1,j+1,k)-f(i-1,j+1,k)-f(i+1,j,k)+f(i-1,j,k))/(2.0*dx*dx);
	pip=0;

	return grad;
}

double gradient::ydxz(fdm* a, field& f)
{
    pip=1;
	grad = (0.5*(f(i+1,j,k+1)+f(i+1,j+1,k+1)) - 0.5*(f(i-1,j,k+1)+f(i-1,j+1,k+1))
	-0.5*(f(i+1,j,k-1)+f(i+1,j+1,k-1)) + 0.5*(f(i-1,j,k-1)+f(i-1,j+1,k-1)))/(4.0*dx*dx);
	pip=0;

	return grad;
}

double gradient::ydyz(fdm* a, field& f)
{
    pip=2;
	grad = (f(i,j+1,k+1)-f(i,j+1,k-1)-f(i,j,k+1)+f(i,j,k-1))/(2.0*dx*dx);
	pip=0;

	return grad;
}

// **********************************************************
// ZX2
// **********************************************************

double gradient::zdx(fdm* a, field& f)
{
    pip=1;
	grad = (0.5*(f(i+1,j,k)+f(i+1,j,k+1)) - 0.5*(f(i-1,j,k)+f(i-1,j,k+1)))/(2.0*dx);
    pip=0;

	return grad;
}

double gradient::zdy(fdm* a, field& f)
{
    pip=2;
	grad = (0.5*(f(i,j+1,k)+f(i,j+1,k+1)) - 0.5*(f(i,j-1,k)+f(i,j-1,k+1)))/(2.0*dx);
    pip=0;

	return grad;
}


double gradient::zdz(fdm* a, field& f)
{
    pip=3;
	grad = (f(i,j,k+1) - f(i,j,k))/(dx);
	pip=0;

	return grad;
}

// **********************************************************
// ZDXX2
// **********************************************************

double gradient::zdxx(fdm* a, field& f)
{
    pip=1;
	grad = (0.5*(f(i-1,j,k)+f(i-1,j,k+1)) - 2.0*0.5*(f(i,j,k)+f(i,j,k+1)) + 0.5*(f(i+1,j,k)+f(i+1,j,k+1)))/(dx*dx);
    pip=0;

	return grad;
}

double gradient::zdyy(fdm* a, field& f)
{
    pip=2;
	grad = (0.5*(f(i,j-1,k)+f(i,j-1,k+1)) - 2.0*0.5*(f(i,j,k)+f(i,j,k+1)) + 0.5*(f(i,j+1,k)+f(i,j+1,k+1)))/(dx*dx);
    pip=0;

	return grad;
}

double gradient::zdzz(fdm* a, field& f)
{
    pip=3;
	grad = (f(i,j,k+2)-f(i,j,k+1)-f(i,j,k)+f(i,j,k-1))/(2.0*dx*dx);
    pip=0;

	return grad;
}

// **********************************************************
// ZDXY2
// **********************************************************

double gradient::zdxy(fdm* a, field& f)
{

	pip=1;
	grad = (0.5*(f(i+1,j+1,k)+f(i+1,j+1,k+1)) - 0.5*(f(i-1,j+1,k)+f(i-1,j+1,k+1))
	-0.5*(f(i+1,j-1,k)+f(i+1,j-1,k+1)) + 0.5*(f(i-1,j-1,k)+f(i-1,j-1,k+1)))/(4.0*dx*dx);
	pip=0;

	return grad;
}

double gradient::zdxz(fdm* a, field& f)
{
    pip=3;
	grad = (f(i+1,j,k+1)-f(i-1,j,k+1)-f(i+1,j,k)+f(i-1,j,k))/(2.0*dx*dx);
	pip=0;

	return grad;
}

double gradient::zdyz(fdm* a, field& f)
{
    pip=3;
	grad = (f(i,j+1,k+1)-f(i,j-1,k+1)-f(i,j+1,k)+f(i,j-1,k))/(2.0*dx*dx);
    pip=0;

	return grad;
}
