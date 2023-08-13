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

#include"fnpf_coastline.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

void fnpf_coastline::disc(lexer *p, ghostcell *pgc, slice &f)
{
    double dx,dy,xmin,xplus,ymin,yplus;
    double lsv,lsSig;
    double dnorm,deltax,sign;
    
    SLICELOOP4
    {
    
    dx=0.0;
	dy=0.0;
    
	lsv=f(i,j);
    lsSig=lsv/sqrt(lsv*lsv);

    if(fabs(lsv)<1.0e-8)
    lsSig=1.0;

// x	
	xmin=(lsv-f(i-1,j))/p->DXP[IM1];
	xplus=(f(i+1,j)-lsv)/p->DXP[IP];
	
	if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
	dx=dswenox(f,1.0);

	if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
	dx=dswenox(f,-1.0);

	if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
	dx=0.0;

// y
    if(p->j_dir==1)
    {
	ymin=(lsv-f(i,j-1))/p->DYP[JM1];
	yplus=(f(i,j+1)-lsv)/p->DYP[JP];
	
	if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
	dy=dswenoy(f,1.0);

	if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
	dy=dswenoy(f,-1.0);

	if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
	dy=0.0;
    }

	dnorm=sqrt(dx*dx + dy*dy);
	
    deltax = 0.5*(p->DXN[IP] + p->DYN[JP]);
    
	sign=lsv/sqrt(lsv*lsv+ dnorm*dnorm*deltax*deltax);

	L(i,j) = -(sign*dnorm - sign);
    }
    
    
}
