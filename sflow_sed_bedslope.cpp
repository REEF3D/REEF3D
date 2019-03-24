
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

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"

void sflow_sediment_f::bedslope(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double uvel,vvel;
    double nx,ny,nz,norm;
    double nx0,ny0;
    double beta;
    
    SLICELOOP4
    {
    // beta
    uvel=0.5*(b->P(i,j)+b->P(i-1,j));

    vvel=0.5*(b->Q(i,j)+b->Q(i,j-1));


	// 1
	if(uvel>0.0 && vvel>0.0 && fabs(uvel)>1.0e-10)
	beta = atan(fabs(vvel/uvel));

	// 2
	if(uvel<0.0 && vvel>0.0 && fabs(vvel)>1.0e-10)
	beta = PI*0.5 + atan(fabs(uvel/vvel));

	// 3
	if(uvel<0.0 && vvel<0.0 && fabs(uvel)>1.0e-10)
	beta = PI + atan(fabs(vvel/uvel));

	// 4
	if(uvel>0.0 && vvel<0.0 && fabs(vvel)>1.0e-10)
	beta = 1.5*PI + atan(fabs(uvel/vvel));

	//------

	if(uvel>0.0 && fabs(vvel)<=1.0e-10)
	beta = 0.0;

	if(fabs(uvel)<=1.0e-10 && vvel>0.0)
	beta = PI*0.5;

	if(uvel<0.0 && fabs(vvel)<=1.0e-10)
	beta = PI;

	if(fabs(uvel)<=1.0e-10 && vvel<0.0)
	beta = PI*1.5;

	if(fabs(uvel)<=1.0e-10 && fabs(vvel)<=1.0e-10)
	beta = 0.0;
    
    
    // n = b x c
    
    //nx = 
   
	/*
    // bed normal
	nx0=(b->depth(i+1,j)-b->depth(i-1,j))/(2.0*p->dx);
    
    if(p->flag4[Im1JK]==OBJ)
    nx0=(b->depth(i+1,j)-b->depth(i,j))/(p->dx);
    
    if(p->flag4[Ip1JK]==OBJ)
    nx0=(b->depth(i,j)-b->depth(i-1,j))/(p->dx);
    
    
	ny0=(b->depth(i,j+1)-b->depth(i,j-1))/(2.0*p->dx);
    
    if(p->flag4[IJm1K]==OBJ)
    ny0=(b->depth(i,j+1)-b->depth(i,j))/(p->dx);
    
    if(p->flag4[IJp1K]==OBJ)
    ny0=(b->depth(i,j)-b->depth(i,j-1))/(p->dx);
    
    
	nz =(b->depth(i,j,k+1)-b->depth(i,j,k-1))/(2.0*p->dx);

	norm=sqrt(nx0*nx0 + ny0*ny0 + nz*nz);
	
	nx0/=norm>1.0e-20?norm:1.0e20;
	ny0/=norm>1.0e-20?norm:1.0e20;
	nz /=norm>1.0e-20?norm:1.0e20;
	
    // rotate bed normal
	
	beta=-beta;
    nx = (cos(beta)*nx0-sin(beta)*ny0);
	ny = (sin(beta)*nx0+cos(beta)*ny0);
    
    teta  = atan(nx/(fabs(nz)>1.0e-15?nz:1.0e20));
    alpha = atan(ny/(fabs(nz)>1.0e-15?nz:1.0e20));

    if(fabs(nx)<1.0e-10 && fabs(ny)<1.0e-10)
    gamma=0.0;

	if(fabs(nx)>=1.0e-10 || fabs(ny)>=1.0e-10)
	gamma = PI*0.5 - acos(	(nx*nx + ny*ny + nz*0.0)/( sqrt(nx*nx + ny*ny + nz*nz )*sqrt(nx*nx + ny*ny + nz*0.0))+1e-20);
	

	phi = midphi + (teta/(fabs(gamma)>1.0e-20?gamma:1.0e20))*delta;*/
    }
}
	