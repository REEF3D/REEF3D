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

#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

void ghostcell::gcsl_start1(lexer *p, slice &f, int gcv)
{
    starttime=timer();
    
	QQGCSL1LOOP
	gcsldistro1(p,f,p->gcbsl1[qq][0], p->gcbsl1[qq][1], p->gcbsl1[qq][5], p->gcdsl1[qq], gcv, p->gcbsl1[qq][4], p->gcbsl1[qq][3]);

	endtime=timer();
	p->gctime+=endtime-starttime;
	
	//  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcslparax(p,f,1);
	gcslparacox(p,f,gcv);
    gcslparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    dgcslpol1(p,f);
}

void ghostcell::gcsl_start2(lexer *p, slice &f, int gcv)
{
    starttime=timer();
    
	QQGCSL2LOOP
	gcsldistro2(p,f,p->gcbsl2[qq][0], p->gcbsl2[qq][1], p->gcbsl2[qq][5], p->gcdsl2[qq], gcv, p->gcbsl2[qq][4], p->gcbsl2[qq][3]);
    
	endtime=timer();
	p->gctime+=endtime-starttime;
	
	//  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcslparax(p,f,2);
	gcslparacox(p,f,gcv);
    gcslparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    dgcslpol2(p,f);
}

void ghostcell::gcsl_start3(lexer *p, slice &f, int gcv)
{
    starttime=timer();
    QQGCSL3LOOP
	gcsldistro3(p,f,p->gcbsl3[qq][0], p->gcbsl3[qq][1], p->gcbsl3[qq][5], p->gcdsl3[qq], gcv, p->gcbsl3[qq][4], p->gcbsl3[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
	
	//  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcslparax(p,f,3);
	gcslparacox(p,f,gcv);
	gcslparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
}

void ghostcell::gcsl_start4(lexer *p, slice &f, int gcv)
{	
	starttime=timer();
	QQGCSL4LOOP
	gcsldistro4(p,f,p->gcbsl4[qq][0],p->gcbsl4[qq][1], p->gcbsl4[qq][5], p->gcdsl4[qq], gcv, p->gcbsl4[qq][4], p->gcbsl4[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
	
	//  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcslparax(p,f,4);
	gcslparacox(p,f,gcv);
	gcslparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    dgcslpol4(p,f);
}

void ghostcell::gcsl_start4int(lexer *p, sliceint &f, int gcv)
{
    starttime=timer();
	QQGCSL4LOOP
	gcsldistro4int(p,f,p->gcbsl4[qq][0],p->gcbsl4[qq][1], p->gcbsl4[qq][5], p->gcdsl4[qq], gcv, p->gcbsl4[qq][4], p->gcbsl4[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
	//  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcslparax_int(p,f,4);
    gcslparacox_int(p,f,gcv);
	gcslparacox_int(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
}

void ghostcell::gcsl_start4a(lexer *p, slice &f, int gcv)
{
    starttime=timer();
	QQGCSL4ALOOP
	gcsldistro4a(p,f,p->gcbsl4a[qq][0], p->gcbsl4a[qq][1], p->gcbsl4a[qq][5], p->gcdsl4a[qq], gcv, p->gcbsl4a[qq][4], p->gcbsl4a[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;
	
	//  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcslparax(p,f,4);
	gcslparacox(p,f,gcv);
    gcslparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
}

void ghostcell::gcsl_start4V(lexer *p, vec2D &f, int gcv, cpt2D &C)
{
	
	starttime=timer();
	QQGCSL4LOOP
	gcsldistro4V(p,f,C,p->gcbsl4[qq][0],p->gcbsl4[qq][1], p->gcbsl4[qq][5], p->gcdsl4[qq], gcv, p->gcbsl4[qq][4], p->gcbsl4[qq][3], p->gcbsl4[qq][5]);
	endtime=timer();
	p->gctime+=endtime-starttime;
    
    //  MPI Boundary Swap
    if(p->M10>0)
    {
    starttime=timer();
	gcparaxvec2D(p,f,4,C);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
	

}
