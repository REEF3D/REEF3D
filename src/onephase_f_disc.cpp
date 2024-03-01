/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"onephase_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"

void onephase_f::disc(lexer *p, fdm* a, field &f)
{	
	/*LOOP
	a->L(i,j,k) = 0.0;
	
    AIRLOOP
	{
	dx=0.0;
	dy=0.0;
	dz=0.0;
	lsv=f(i,j,k);
    lsSig=lsv/sqrt(lsv*lsv);

    if(fabs(lsv)<1.0e-8)
    lsSig=1.0;

// x
	xmin=(lsv-f(i-1,j,k))/p->DXP[IM1];
	xplus=(f(i+1,j,k)-lsv)/p->DXP[IP];
	
	if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
	dx=ddwenox(a,f,1.0);

	if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
	dx=ddwenox(a,f,-1.0);

	if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
	dx=0.0;

// y
	ymin=(lsv-f(i,j-1,k))/p->DYP[JM1];
	yplus=(f(i,j+1,k)-lsv)/p->DYP[JP];
	
	if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
	dy=ddwenoy(a,f,1.0);

	if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
	dy=ddwenoy(a,f,-1.0);

	if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
	dy=0.0;

// z
	zmin=(lsv-f(i,j,k-1))/p->DZP[KM1];
	zplus=(f(i,j,k+1)-lsv)/p->DZP[KP];
	
	if(zmin*lsSig>0.0 && zplus*lsSig>-zmin*lsSig)
	dz=ddwenoz(a,f,1.0);

	if(zplus*lsSig<0.0 && zmin*lsSig<-zplus*lsSig)
	dz=ddwenoz(a,f,-1.0);

	if(zplus*lsSig>0.0 && zmin*lsSig<0.0)
	dz=0.0;
					

	dnorm=sqrt(dx*dx + dy*dy + dz*dz);
    
    if(p->j_dir==0)
    deltax = (1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
	
    if(p->j_dir==1)
    deltax = (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
	
	sign=lsv/sqrt(lsv*lsv+ dnorm*dnorm*deltax*deltax);

	a->L(i,j,k) = -(sign*dnorm - sign);
    }*/
}