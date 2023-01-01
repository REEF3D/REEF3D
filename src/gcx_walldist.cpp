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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"reini.h"
#include"reini_walld.h"


void ghostcell::walldistance(lexer *p, fdm *a, ghostcell *pgc, convection *pdisc, reini *preini, ioflow *pflow,  field& walldist)
{
	int ic,jc,kc;
	double xdist,ydist,zdist;
	
	MALOOP
    walldist(i,j,k)=0.0;
	
	LOOP
    walldist(i,j,k)=1.0e9;
	
	pgc->gcparax(p,walldist,4);
	
	GC4LOOP
	{
	ic=p->gcb4[n][0];
	jc=p->gcb4[n][1];
	kc=p->gcb4[n][2];	

	
	if(p->gcb4[n][4]==5 || p->gcb4[n][4]==21 || p->gcb4[n][4]==22)
	{
		
		if(p->gcb4[n][3]==1)
		for(i=ic;i<p->knox;++i)
		{
		j=jc;
		k=kc;
		xdist = fabs(double(i-ic)*p->DXM) + 0.5*p->DXM;
		PCHECK
		walldist(i,j,k)=MIN(walldist(i,j,k),xdist);
		}
		
		if(p->gcb4[n][3]==3)
		for(j=jc;j<p->knoy;++j)
		{
		i=ic;
		k=kc;
		ydist = fabs(double(j-jc)*p->DXM) + 0.5*p->DXM;
		PCHECK
		walldist(i,j,k)=MIN(walldist(i,j,k),ydist);
		}
		
		if(p->gcb4[n][3]==5)
		for(k=kc;k<p->knoz;++k)
		{
		i=ic;
		j=jc;
		zdist = fabs(double(k-kc)*p->DXM) + 0.5*p->DXM;
		PCHECK
		walldist(i,j,k)=MIN(walldist(i,j,k),zdist);
		}
		
		if(p->gcb4[n][3]==4)
		for(i=0;i<=ic;++i)
		{
		j=jc;
		k=kc;
		xdist = fabs(double(i-ic)*p->DXM) + 0.5*p->DXM;
		PCHECK
		walldist(i,j,k)=MIN(walldist(i,j,k),xdist);
		}
		
		if(p->gcb4[n][3]==2)
		for(j=0;j<=jc;++j)
		{
		i=ic;
		k=kc;
		ydist = fabs(double(j-jc)*p->DXM) + 0.5*p->DXM;
		PCHECK
		walldist(i,j,k)=MIN(walldist(i,j,k),ydist);
		}
		
		if(p->gcb4[n][3]==6)
		for(k=0;k<=kc;++k)
		{
		i=ic;
		j=jc;
		zdist = fabs(double(k-kc)*p->DXM) + 0.5*p->DXM;
		PCHECK
		walldist(i,j,k)=MIN(walldist(i,j,k),zdist);
		}
	}
	}
	//LOOP
	//walldist(i,j,k)=1.0;
	
	pgc->gcparax(p,walldist,4);
	gcparacox(p,walldist,4);
	gcparacox(p,walldist,4);
	
	
	// calculate global position of gcb cell
	count=0;

    GC4LOOP
    if(p->gcb4[n][4]==5|| p->gcb4[n][4]==21 || p->gcb4[n][4]==22)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	
	if(p->gcb4[n][3]==1)
	walldist(i-1,j,k)=-0.5*p->DXM;  

	if(p->gcb4[n][3]==4)
	walldist(i+1,j,k)=-0.5*p->DXM;   

	if(p->gcb4[n][3]==3)
	walldist(i,j-1,k)=-0.5*p->DXM;   
	
	if(p->gcb4[n][3]==2)
	walldist(i,j+1,k)=-0.5*p->DXM; 
	
	if(p->gcb4[n][3]==5)
	walldist(i,j,k-1)=-0.5*p->DXM; 
	
	if(p->gcb4[n][3]==6)
	walldist(i,j,k+1)=-0.5*p->DXM; 
    }
	
	
	reini_walld reini(p,a);
	
	reini.start(a,p,walldist,pgc,pflow);

	pgc->gcparax(p,walldist,4);
	gcparacox(p,walldist,4);
	gcparacox(p,walldist,4);
	
	GC4LOOP
    if(p->gcb4[n][4]==5|| p->gcb4[n][4]==21 || p->gcb4[n][4]==22)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	
	if(p->gcb4[n][3]==1)
	walldist(i-1,j,k)=0.0;  

	if(p->gcb4[n][3]==4)
	walldist(i+1,j,k)=0.0;   

	if(p->gcb4[n][3]==3)
	walldist(i,j-1,k)=0.0;   
	
	if(p->gcb4[n][3]==2)
	walldist(i,j+1,k)=0.0; 
	
	if(p->gcb4[n][3]==5)
	walldist(i,j,k-1)=0.0; 
	
	if(p->gcb4[n][3]==6)
	walldist(i,j,k+1)=0.0; 
    }
	
}


