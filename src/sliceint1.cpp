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

#include"sliceint1.h"
#include"lexer.h"
#include"fdm.h"

sliceint1::sliceint1(lexer *p)
{
    imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    
	fieldalloc(p);
	fieldgcalloc(p);
	
	pp=p;
}

sliceint1::~sliceint1()
{
    int a,b;
    
    for(a=0;a<gcfeldsize;++a)
    for(b=0;b<4;++b)
	delete [ ] gcfeld[a][b];
    
	for(a=0;a<gcfeldsize;++a)
	delete [ ] gcfeld[a];

	delete [ ] gcfeld;

	delete [ ] V;
}

void sliceint1::fieldalloc(lexer* p)
{
	int gridsize = imax*jmax;
	p->Iarray(V,gridsize);
}

void sliceint1::dealloc(lexer* p)
{
    int a,b;
    
    for(a=0;a<gcfeldsize;++a)
    for(b=0;b<4;++b)
	delete [ ] gcfeld[a][b];
    
	for(a=0;a<gcfeldsize;++a)
	delete [ ] gcfeld[a];

	delete [ ] gcfeld;

	delete [ ] V;
    
    gcfeldsize=feldsize=0;
}

void sliceint1::resize(lexer* p)
{
    if(p->gcsl_extra1*p->margin>gcfeldsize)
    cout<<p->mpirank<<" Slice4 Resize: "<<gcfeldsize<<" "<<p->gcsl_extra1*p->margin<<endl;
}

void sliceint1::fieldgcalloc(lexer* p)
{
    gcfeldsize=p->gcsl_extra1*p->margin;
	gcsl_extra=gcfeldsize;
    
	gcfeldsize+=(p->gcbsl1_count);
	
	p->Iarray(gcfeld,gcfeldsize,4,4);
}

int & sliceint1::operator()(int ii, int jj)
{			
	if(pp->mgcsl1[(ii-imin)*jmax + (jj-jmin)]<2)
	return V[(ii-imin)*jmax + (jj-jmin)];
	
	
	iter=(ii-imin)*jmax + (jj-jmin);
	
		di=ii-i;
		dj=jj-j;

		if(pip==4)
		return V[iter];
		
		if(di==0 && dj==0 && pip!=1)
		return V[iter];
        
        if(di==0 && dj==0 && pip==1)
		di=1;

	  
//1
		if(di<0 && (dj==0||pip==1))
		{
			if(pp->gcslorig1[pp->mgcsl1[iter]-10][0][-di]==0)
            {
            if(di<-2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][0][-di-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][0][-di-1];
            
            if(di<-2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][0][-di-2]==1)
			return gcfeld[pp->mgcsl1[iter]-10][0][-di-2];
            
            if(di<-1)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][0][-di-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][0][-di-1];
            
            return V[iter];
            }
			
			if(pp->gcslorig1[pp->mgcsl1[iter]-10][0][-di]==1)
			return gcfeld[pp->mgcsl1[iter]-10][0][-di];
		}
//4
		if(di>0 && (dj==0||pip==1))
		{
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][3][di]==0)
            {
            if(di>2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][3][di-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][3][di-1];
            
            if(di>2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][3][di-2]==1)
			return gcfeld[pp->mgcsl1[iter]-10][3][di-2];
            
            if(di>1)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][3][di-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][3][di-1];
            
            return V[iter];
            }
			
			if(pp->gcslorig1[pp->mgcsl1[iter]-10][3][di]==1)
			return gcfeld[pp->mgcsl1[iter]-10][3][di];
		}

//3
		if(dj<0 && (di==0||pip==2))
		{
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][2][-dj]==0)
            {
            if(dj<-2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][2][-dj-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][2][-dj-1];
            
            if(dj<-2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][2][-dj-2]==1) 
			return gcfeld[pp->mgcsl1[iter]-10][2][-dj-2];
            
            if(dj<-1)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][2][-dj-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][2][-dj-1];
            
            return V[iter];
            }
			
			if(pp->gcslorig1[pp->mgcsl1[iter]-10][2][-dj]==1)
			return gcfeld[pp->mgcsl1[iter]-10][2][-dj];
		}
//2
		if(dj>0 && (di==0||pip==2))
		{
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][1][dj]==0)
            {
            if(dj>2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][1][dj-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][1][dj-1];
            
            if(dj>2)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][1][dj-2]==1) 
			return gcfeld[pp->mgcsl1[iter]-10][1][dj-2];
            
            if(dj>1)
            if(pp->gcslorig1[pp->mgcsl1[iter]-10][1][dj-1]==1)
			return gcfeld[pp->mgcsl1[iter]-10][1][dj-1];
            
            return V[iter];
            }
			
			if(pp->gcslorig1[pp->mgcsl1[iter]-10][1][dj]==1)
			return gcfeld[pp->mgcsl1[iter]-10][1][dj];
		}


	return V[iter];
	
}

