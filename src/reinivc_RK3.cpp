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

#include"reinivc_RK3.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"

reinivc_RK3::reinivc_RK3(lexer* p):gradient(p),iflag(p),
												epsi(p->F45*p->DXM)
{
	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

    gcval_iniphi=50;

    if((p->F61>1.0e-20 || p->F60>1.0e-20) && p->F50==1)
    gcval_iniphi=51;
	
	if((p->F62>1.0e-20 || p->F60>1.0e-20) && p->F50==2)
    gcval_iniphi=52;

	gcval_ro=1;
	dt=p->F43*p->DXM;

	if(p->F46==1)
	ppicard = new picard_f(p);

	if(p->F46!=1)
	ppicard = new picard_void(p);
}

reinivc_RK3::~reinivc_RK3()
{
}

void reinivc_RK3::start(fdm* a,lexer* p,field& b, ghostcell* pgc,ioflow* pflow)
{ 
    field4 brk1(p),brk2(p);
    
    starttime=pgc->timer();
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<"initializing level set..."<<endl<<endl;
	reiniter=2*int(p->maxlength/(dt));
	pgc->start4(p,b,gcval_iniphi);
	}

	if(p->count>0)
	step(a,p);

	pflow->fsfrkin(p,a,pgc,brk1);
    pflow->fsfrkin(p,a,pgc,brk2);
    pflow->fsfrkout(p,a,pgc,brk1);
    pflow->fsfrkout(p,a,pgc,brk2);


    ppicard->volcalc(p,a,pgc,a->phi);
    pgc->start4(p,b,gcval_phi);
	
	if(p->count>0)
	interface_cells(p,a,b);

    for(int q=0;q<reiniter;++q)
    {
	// Step 1
	disc(p,a,b);

	LOOP
	brk1(i,j,k)=b(i,j,k)+ dt*a->L(i,j,k);

	if(p->count==0)
	pgc->start4(p,brk1,gcval_iniphi);

	pgc->start4(p,brk1,gcval_phi);

    // Step 2
    disc(p,a,brk1);

	LOOP
	brk2(i,j,k)= 0.75*b(i,j,k) + 0.25*brk1(i,j,k) + 0.25*dt*a->L(i,j,k);

	if(p->count==0)
	pgc->start4(p,brk2,gcval_iniphi);

	pgc->start4(p,brk2,gcval_phi);

    // Step 3
    disc(p,a,brk2);

	LOOP
	b(i,j,k)=(1.0/3.0)*b(i,j,k) + (2.0/3.0)*brk2(i,j,k) + (2.0/3.0)*dt*a->L(i,j,k);

	if(p->count==0)
	pgc->start4(p,b,gcval_iniphi);

	pgc->start4(p,b,gcval_phi);
	}
	
	if(p->count>0)
	correction(p,a,b);
	
	pgc->start4(p,b,gcval_phi);

	ppicard->correct_ls(p,a,pgc,a->phi);
	
	if(p->count>0)
	finalize(p,a);

	p->reinitime=pgc->timer()-starttime;
}

void reinivc_RK3::startV(fdm* a,lexer* p,vec &f, ghostcell* pgc,ioflow* pflow)
{ 
    
}

void reinivc_RK3::interface_cells(lexer *p, fdm* a, field& b)
{
	int check;
	double zero=0.0;
	
	count=0;
	LOOP
	iflag(i,j,k)=0;
	
	LOOP
    {
        check=1;

        if(b(i,j,k)<zero && b(i+1,j,k)<zero && b(i+1,j+1,k)<zero && b(i,j+1,k)<zero &&
           b(i,j,k+1)<zero && b(i+1,j,k+1)<zero && b(i+1,j+1,k+1)<zero && b(i,j+1,k+1)<zero)
        check=0;
		
		if(b(i,j,k)>zero && b(i+1,j,k)>zero && b(i+1,j+1,k)>zero && b(i,j+1,k)>zero &&
           b(i,j,k+1)>zero && b(i+1,j,k+1)>zero && b(i+1,j+1,k+1)>zero && b(i,j+1,k+1)>zero)
        check=0;
		
		if(check==1)
		{
		iflag(i,j,k)=1;
		++count;
		}
	}
	
	numvert=count;
	
	p->Iarray(ijk,numvert,3);
    p->Darray(ls,numvert);
	p->Darray(ls0,numvert);
	p->Darray(ls1,numvert);
	
	count=0;
	LOOP
    {
		if(iflag(i,j,k)==1)
		{
		ijk[count][0] = i;
		ijk[count][1] = j;
		ijk[count][2] = k;	
		ls[count] = b(i,j,k);
		ls0[count] = b(i,j,k);
				
		++count;
		}
	}
}

void reinivc_RK3::correction(lexer *p, fdm* a, field& b)
{

	double phival,dval,H0,denom;

	for(n=0;n<numvert;++n)
	ls1[n]=ls[n];
	
	dV=pow(p->DXM,3.0);
	
	for(n=0;n<numvert;++n)
    {
		dV1 = 1.0;
		eta = 0.0;
		
		count=0;
		while(dV1>1.0e-10*dV && count<1000)
		{
			
			phival = ls0[n];
			if(phival>epsi)
			H0=1.0;

			if(phival<-epsi)
			H0=0.0;

			if(fabs(phival)<=epsi)
			H0=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
			
			dval = ls[n];
			if(dval>epsi)
			H=1.0;

			if(dval<-epsi)
			H=0.0;

			if(fabs(dval)<=epsi)
			H=0.5*(1.0 + dval/epsi + (1.0/PI)*sin((PI*dval)/epsi));

			
			dV1 = dV*(H0-H);
			
			denom = fabs(ls[n])>1.0e-19?fabs(ls[n]):1.0e20;
			eta = dV1/denom;
			
			ls[n] += eta;
			
			
			++count;
		}
	}
	
	
	for(n=0;n<numvert;++n)
    {
    i=ijk[n][0];
    j=ijk[n][1];
    k=ijk[n][2];
    b(i,j,k)=ls[n];
    }

}

void reinivc_RK3::disc(lexer *p, fdm* a, field& b)
{	
	LOOP
	a->L(i,j,k) = 0.0;
	
    LOOP
	{
	dx=0.0;
	dy=0.0;
	dz=0.0;
	lsv=b(i,j,k);
    lsSig=lsv/sqrt(lsv*lsv);

    if(fabs(lsv)<1.0e-8)
    lsSig=1.0;

// x
	xmin=(lsv-b(i-1,j,k))/p->DXP[IM1];
	xplus=(b(i+1,j,k)-lsv)/p->DXP[IP];
	
	if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
	dx=ddwenox(a,b,1.0);

	if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
	dx=ddwenox(a,b,-1.0);

	if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
	dx=0.0;

// y
	ymin=(lsv-b(i,j-1,k))/p->DYP[JM1];
	yplus=(b(i,j+1,k)-lsv)/p->DYP[JP];
	
	if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
	dy=ddwenoy(a,b,1.0);

	if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
	dy=ddwenoy(a,b,-1.0);

	if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
	dy=0.0;

// z
	zmin=(lsv-b(i,j,k-1))/p->DZP[KM1];
	zplus=(b(i,j,k+1)-lsv)/p->DZP[KP];
	
	if(zmin*lsSig>0.0 && zplus*lsSig>-zmin*lsSig)
	dz=ddwenoz(a,b,1.0);

	if(zplus*lsSig<0.0 && zmin*lsSig<-zplus*lsSig)
	dz=ddwenoz(a,b,-1.0);

	if(zplus*lsSig>0.0 && zmin*lsSig<0.0)
	dz=0.0;
					

	dnorm=sqrt(dx*dx + dy*dy + dz*dz);
    
    if(p->j_dir==0)
    deltax = (1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
	
    if(p->j_dir==1)
    deltax = (1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
	
	sign=lsv/sqrt(lsv*lsv+ dnorm*dnorm*deltax*deltax);

	a->L(i,j,k) = -(sign*dnorm - sign);
    }
}

void reinivc_RK3::step(fdm* a, lexer* p)
{
	dt=p->F43*deltax;

	reiniter=p->F44;
}

void reinivc_RK3::finalize(lexer* p, fdm* a)
{
	p->del_Iarray(ijk,numvert,3);
    p->del_Darray(ls,numvert);
	p->del_Darray(ls0,numvert);
	p->del_Darray(ls1,numvert);
}
