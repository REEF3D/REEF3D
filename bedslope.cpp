/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"bedslope.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fnpf_weno.h"

bedslope::bedslope(lexer *p) : norm_vec(p)
{
    midphi=p->S81*(PI/180.0);
    delta=p->S82*(PI/180.0);
    
     pdx = new fnpf_weno(p);
}

bedslope::~bedslope()
{
}

void bedslope::slope(lexer *p, fdm * a, ghostcell *pgc, double &teta, double &alpha, double &gamma, double &phi)
{
    double uvel,vvel;
    double nx,ny,nz,norm;
    double nx0,ny0;
    double nz0,bx0,by0;
    
    // beta
    uvel=0.5*(a->P(i,j)+a->P(i-1,j));

    vvel=0.5*(a->Q(i,j)+a->Q(i,j-1));


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
   
    
    // ----
    
     bx0 = (a->bedzh(i+1,j)-a->bedzh(i-1,j))/(2.0*p->DXM);
     by0 = (a->bedzh(i,j+1)-a->bedzh(i,j-1))/(2.0*p->DXM);
     /*
     signx=signy=1.0;
     
     if(bx0<0.0)
     signx=-1.0;
     
     if(by0<0.0)
     signy=-1.0;
     
     bx0 = pdx->sx(p,a->bedzh,signx);
     bx0 = pdx->sy(p,a->bedzh,signy);
     */
     nx0 = bx0/sqrt(bx0*bx0 + by0*by0 + 1.0);
     ny0 = by0/sqrt(bx0*bx0 + by0*by0 + 1.0);
     nz0 = 1.0;
     
     norm=sqrt(nx0*nx0 + ny0*ny0 + nz0*nz0);
     
     
    nx0/=norm>1.0e-20?norm:1.0e20;
	ny0/=norm>1.0e-20?norm:1.0e20;
	nz0/=norm>1.0e-20?norm:1.0e20;
	
    // rotate bed normal
	beta=-beta;
    nx = (cos(beta)*nx0-sin(beta)*ny0);
	ny = (sin(beta)*nx0+cos(beta)*ny0);
    nz = nz0;
    
    teta  = -atan(nx/(fabs(nz)>1.0e-15?nz:1.0e20));
    alpha =  fabs(atan(ny/(fabs(nz)>1.0e-15?nz:1.0e20)));
    
    
    //-----------

    if(fabs(nx)<1.0e-10 && fabs(ny)<1.0e-10)
    gamma=0.0;

	if(fabs(nx)>=1.0e-10 || fabs(ny)>=1.0e-10)
	gamma = PI*0.5 - acos(	(nx*nx + ny*ny + nz*0.0)/( sqrt(nx*nx + ny*ny + nz*nz )*sqrt(nx*nx + ny*ny + nz*0.0))+1e-20);
	
    
    phi = midphi + (teta/(fabs(gamma)>1.0e-20?fabs(gamma):1.0e20))*delta; 

}

