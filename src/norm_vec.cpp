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

#include"norm_vec.h"
#include"lexer.h"
#include"fdm.h"

norm_vec::norm_vec(lexer* p) : ddweno_f_nug(p), dx(p->DXM)
{
}

norm_vec::~norm_vec()
{
}

double norm_vec::normvec_x(fdm *a, field &f)
{
    double lsv,lsSig;
    double xmin;
    double xplus;
    double nx;

	nx=0.0;

	lsv=f(i,j,k);
	lsSig=lsv/sqrt(lsv*lsv+dx*dx);

	xmin=(lsv-f(i-1,j,k))/dx;
	xplus=(f(i+1,j,k)-lsv)/dx;


// x
	if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
	nx=ddwenox(f,1.0);

	if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
	nx=ddwenox(f,-1.0);

	if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
	nx=0.5*ddwenox(f,1.0) + 0.5*ddwenox(f,-1.0);

	nx=fabs(nx)>1.0e-15?nx:1.0e-15;

    return nx;

}

double norm_vec::normvec_y(fdm *a, field &f)
{
    double lsv,lsSig;
    double ymin;
    double yplus;
    double ny;

	ny=0.0;

	lsv=f(i,j,k);
	lsSig=lsv/sqrt(lsv*lsv+dx*dx);

	ymin=(lsv-f(i,j-1,k))/dx;
	yplus=(f(i,j+1,k)-lsv)/dx;

// y
	if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
	ny=ddwenoy(f,1.0);

	if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
	ny=ddwenoy(f,-1.0);

	if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
	ny=0.5*ddwenoy(f,1.0) + 0.5*ddwenoy(f,-1.0);

	ny=fabs(ny)>1.0e-15?ny:1.0e-15;

    return ny;
}

double norm_vec::normvec_z(fdm *a, field &f)
{
    double lsv,lsSig;
    double zmin;
    double zplus;
    double nz;

	nz=0.0;
	lsv=f(i,j,k);
	lsSig=lsv/sqrt(lsv*lsv+dx*dx);

	zmin=(lsv-f(i,j,k-1))/dx;
	zplus=(f(i,j,k+1)-lsv)/dx;

// z
	if(zmin*lsSig>0.0 && zplus*lsSig>-zmin*lsSig)
	nz=ddwenoz(f,1.0);

	if(zplus*lsSig<0.0 && zmin*lsSig<-zplus*lsSig)
	nz=ddwenoz(f,-1.0);

	if(zplus*lsSig>0.0 && zmin*lsSig<0.0)
	nz=0.5*ddwenoz(f,1.0) + 0.5*ddwenoz(f,-1.0);

	nz=fabs(nz)>1.0e-15?nz:1.0e-15;

	return nz;

}

