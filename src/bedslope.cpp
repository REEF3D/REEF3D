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

#include"bedslope.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"ddweno_f_nug.h"

bedslope::bedslope(lexer *p) : norm_vec(p)
{
    midphi=p->S81*(PI/180.0);
    delta=p->S82*(PI/180.0);
    
    pdx = new ddweno_f_nug(p);
}

bedslope::~bedslope()
{
}

void bedslope::slope_weno(lexer *p, ghostcell *pgc, sediment_fdm *s, field &topo)
{
    double uvel,vvel;
    double nx,ny,nz,norm;
    double nx0,ny0,nx1,ny1,nz1;
    double nz0,bx0,by0;
    
    SLICELOOP4
    {
    
    // beta
    uvel=0.5*(s->P(i,j)+s->P(i-1,j));
    vvel=0.5*(s->Q(i,j)+s->Q(i,j-1));

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
    
    // bed normal
	nx0=-(topo(i+1,j,k)-topo(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);
    
    if(p->flag4[Im1JK]<=SOLID)
    nx0=-(topo(i+1,j,k)-topo(i,j,k))/(p->DXP[IP]);
    
    if(p->flag4[Ip1JK]<=SOLID)
    nx0=-(topo(i,j,k)-topo(i-1,j,k))/(p->DXP[IM1]);
    
    
	ny0=-(topo(i,j+1,k)-topo(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]);
    
    if(p->flag4[IJm1K]<=SOLID)
    ny0=-(topo(i,j+1,k)-topo(i,j,k))/(p->DYP[JP]);
    
    if(p->flag4[IJp1K]<=SOLID)
    ny0=-(topo(i,j,k)-topo(i,j-1,k))/(p->DYP[JM1]);
    
    
	nz0 = (topo(i,j,k+1)-topo(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
    
    nx1 = -pdx->ddwenox(topo, nx0);
    ny1 = -pdx->ddwenoy(topo, ny0);
    nz1 =  pdx->ddwenoz(topo, nz0);
    

	norm=sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1);
	
	nx1/=norm>1.0e-20?norm:1.0e20;
	ny1/=norm>1.0e-20?norm:1.0e20;
	nz1/=norm>1.0e-20?norm:1.0e20;
   
    
    // rotate bed normal
	beta=-beta;
    nx = (cos(beta)*nx1-sin(beta)*ny1);
	ny = (sin(beta)*nx1+cos(beta)*ny1);
    nz = nz1;
    
    s->beta(i,j) = -beta;
    
    s->teta(i,j)  = -atan(nx/(fabs(nz)>1.0e-15?nz:1.0e20));
    s->alpha(i,j) =  fabs(atan(ny/(fabs(nz)>1.0e-15?nz:1.0e20)));
    
    
    //-----------

    if(fabs(nx)<1.0e-10 && fabs(ny)<1.0e-10)
    s->gamma(i,j)=0.0;

	if(fabs(nx)>=1.0e-10 || fabs(ny)>=1.0e-10)
	s->gamma(i,j) = PI*0.5 - acos(	(nx1*nx1 + ny1*ny1 + nz1*0.0)/( sqrt(nx1*nx1 + ny1*ny1 + nz1*nz1 )*sqrt(nx1*nx1 + ny1*ny1 + nz1*0.0))+1e-20);
	
    
    s->phi(i,j) = midphi + MIN(1.0,fabs(s->teta(i,j)/midphi))*(s->teta(i,j)/(fabs(s->gamma(i,j))>1.0e-20?fabs(s->gamma(i,j)):1.0e20))*delta; 
    }
}

void bedslope::slope_cds(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    double uvel,vvel;
    double nx,ny,nz,norm;
    double nx0,ny0;
    double nz0,bx0,by0;
    
    SLICELOOP4
    {
    
    // beta
    uvel=0.5*(s->P(i,j)+s->P(i-1,j));

    vvel=0.5*(s->Q(i,j)+s->Q(i,j-1));


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
    
     bx0 = (s->bedzh(i+1,j)-s->bedzh(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
     by0 = (s->bedzh(i,j+1)-s->bedzh(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
     
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
  
    s->beta(i,j) = -beta;
    
    s->teta(i,j)  = -atan(nx/(fabs(nz)>1.0e-15?nz:1.0e20));
    s->alpha(i,j) =  fabs(atan(ny/(fabs(nz)>1.0e-15?nz:1.0e20)));
    
    //-----------

    if(fabs(nx)<1.0e-10 && fabs(ny)<1.0e-10)
    s->gamma(i,j)=0.0;

    s->gamma(i,j) = atan(sqrt(bx0*bx0 + by0*by0));

    s->phi(i,j) = midphi + MIN(1.0,fabs(s->teta(i,j)/midphi))*(s->teta(i,j)/(fabs(s->gamma(i,j))>1.0e-20?fabs(s->gamma(i,j)):1.0e20))*delta; 
    }
}
