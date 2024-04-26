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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"slice.h"

void sixdof_obj::reini_2D(lexer *p, ghostcell *pgc, slice &f)
{
	if(p->count==0)
	{
	reiniter=2*int(p->maxlength/(1.0*p->DXM));
    
    if(p->mpirank==0)
	cout<<"initializing 6DOF... "<<endl<<endl;
	}

	if(p->count>0)
	reiniter=p->S37;

    for(int q=0;q<reiniter;++q)
    {
	// Step 1
    disc_2D(p,pgc,f);
    
	SLICELOOP4
	lrk1(i,j) = f(i,j) + dts(i,j)*K(i,j);

	pgc->gcsl_start4(p,lrk1,50);
    

    // Step 2
    disc_2D(p,pgc,lrk1);
    
	SLICELOOP4
	lrk2(i,j)=  0.75*f(i,j) + 0.25*lrk1(i,j) + 0.25*dts(i,j)*K(i,j);

	pgc->gcsl_start4(p,lrk2,50);


    // Step 3
    disc_2D(p,pgc,lrk2);
    
	SLICELOOP4
	f(i,j) = (1.0/3.0)*f(i,j) + (2.0/3.0)*lrk2(i,j) + (2.0/3.0)*dts(i,j)*K(i,j);

	pgc->gcsl_start4(p,f,50);
	}
    
    pgc->gcsl_start4(p,f,50);
}


void sixdof_obj::disc_2D(lexer *p, ghostcell *pgc, slice &f)
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

	K(i,j) = -(sign*dnorm - sign);
    }
}


void sixdof_obj::time_preproc_2D(lexer* p)
{	
    n=0;
	SLICELOOP4
	{
	dts(i,j) = p->F43*MIN(p->DXP[IP],p->DYP[JP]);
	++n;
	}
}

