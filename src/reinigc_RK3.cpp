/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
#include"reinigc_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"


reinigc_RK3::reinigc_RK3(lexer *p, fdm *a):gradient(p),deltax(p->DXM),d0(p),wallf(p),epsi(p->F45*p->DXM)
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

    if(p->F61>1.0e-20 || p->F62>1.0e-20)
    gcval_iniphi=gcval_phi;

	gcval_ro=1;
	dt=p->F43*deltax;

	if(p->F46==1)
	ppicard = new picard_f(p);

	if(p->F46!=1)
	ppicard = new picard_void(p);
	
	dV = p->DXM*p->DXM*p->DXM;
}

reinigc_RK3::~reinigc_RK3()
{
}

void reinigc_RK3::start(fdm* a,lexer* p,field& b, ghostcell* pgc,ioflow* pflow)
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
	
	FLUIDLOOP
	d0(i,j,k)=b(i,j,k);
	pgc->start4(p,d0,gcval_phi);
	

    for(int q=0;q<reiniter;++q)
    {
	
	// Step 1
	disc(p,a,pgc,b);

	FLUIDLOOP
	brk1(i,j,k)=b(i,j,k)+ dt*a->L(i,j,k);

	if(p->count==0)
	pgc->start4(p,brk1,gcval_iniphi);

	pgc->start4(p,brk1,gcval_phi);
	
    // Step 2
    disc(p,a,pgc,brk1);

	FLUIDLOOP
	brk2(i,j,k)= 0.75*b(i,j,k) + 0.25*brk1(i,j,k) + 0.25*dt*a->L(i,j,k);

	if(p->count==0)
	pgc->start4(p,brk2,gcval_iniphi);

	pgc->start4(p,brk2,gcval_phi);

    // Step 3
    disc(p,a,pgc,brk2);

	FLUIDLOOP
	b(i,j,k)=(1.0/3.0)*b(i,j,k) + (2.0/3.0)*brk2(i,j,k) + (2.0/3.0)*dt*a->L(i,j,k);

	if(p->count==0)
	pgc->start4(p,b,gcval_iniphi);

	pgc->start4(p,b,gcval_phi);

	}
	
	/*if(p->count>0)
	constraint(p,a,pgc,b);
	pgc->start4(p,b,gcval_phi);*/
	
	ppicard->correct_ls(p,a,pgc,a->phi);

	p->reinitime=pgc->timer()-starttime;
}

void reinigc_RK3::startV(fdm* a,lexer* p,vec &f, ghostcell* pgc,ioflow* pflow)
{ 
    
}

void reinigc_RK3::disc(lexer *p, fdm* a, ghostcell *pgc, field& b)
{	
	FLUIDLOOP
	a->L(i,j,k) = 0.0;
	
    FLUIDLOOP
	if((b(i,j,k)>=0.0 && b(i+1,j,k)>=0.0 && b(i-1,j,k)>=0.0 && b(i,j+1,k)>=0.0 && b(i,j-1,k)>=0.0 && b(i,j,k+1)>=0.0 && b(i,j,k-1)>=0.0) 
	|| (b(i,j,k)<0.0 && b(i+1,j,k)<0.0 && b(i-1,j,k)<0.0 && b(i,j+1,k)<0.0 && b(i,j-1,k)<0.0 && b(i,j,k-1)<0.0 && b(i,j,k+1)<0.0) 
	|| p->count==0 || p->F49==1)
	{
	dstx=0.0;
	dsty=0.0;
	dstz=0.0;
	lsv=b(i,j,k);
    lsSig=lsv/sqrt(lsv*lsv);

    if(fabs(lsv)<1.0e-8)
    lsSig=1.0;

	
// x
	xmin=(lsv-b(i-1,j,k))/deltax;
	xplus=(b(i+1,j,k)-lsv)/deltax;
	
	if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
	dstx=ddwenox(a,b,1.0);

	if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
	dstx=ddwenox(a,b,-1.0);

	if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
	dstx=0.0;

// y
	ymin=(lsv-b(i,j-1,k))/deltax;
	yplus=(b(i,j+1,k)-lsv)/deltax;
	
	if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
	dsty=ddwenoy(a,b,1.0);

	if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
	dsty=ddwenoy(a,b,-1.0);

	if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
	dsty=0.0;

// z
	zmin=(lsv-b(i,j,k-1))/deltax;
	zplus=(b(i,j,k+1)-lsv)/deltax;
	
	if(zmin*lsSig>0.0 && zplus*lsSig>-zmin*lsSig)
	dstz=ddwenoz(a,b,1.0);

	if(zplus*lsSig<0.0 && zmin*lsSig<-zplus*lsSig)
	dstz=ddwenoz(a,b,-1.0);

	if(zplus*lsSig>0.0 && zmin*lsSig<0.0)
	dstz=0.0;


	dnorm=sqrt(dstx*dstx + dsty*dsty + dstz*dstz);

	sign=lsv/sqrt(lsv*lsv+ dnorm*dnorm*deltax*deltax);

	a->L(i,j,k) = -(sign*dnorm - sign);
    }
}

void reinigc_RK3::constraint(lexer *p, fdm* a, ghostcell *pgc, field& b)
{
	
	
	FLUIDLOOP
	{
  
        dx=0.0;
        dy=0.0;
        dz=0.0;
        lsv=b(i,j,k);
        lsSig=lsv/sqrt(lsv*lsv);

        if(fabs(lsv)<1.0e-8)
        lsSig=1.0;
        
        dirac=0.0;
		
		if(fabs(b(i,j,k))<epsi)
		dirac = (1.0/(2.0*epsi)) + (1.0/(2.0*epsi))*cos((PI*b(i,j,k))/epsi);

	
    // x
        xmin=(lsv-b(i-1,j,k))/deltax;
        xplus=(b(i+1,j,k)-lsv)/deltax;
        
        if(xmin*lsSig>0.0 && xplus*lsSig>-xmin*lsSig)
        dx=ddwenox(a,b,1.0);

        if(xplus*lsSig<0.0 && xmin*lsSig<-xplus*lsSig)
        dx=ddwenox(a,b,-1.0);

        if(xplus*lsSig>0.0 && xmin*lsSig<0.0)
        dx=0.0;

    // y
        ymin=(lsv-b(i,j-1,k))/deltax;
        yplus=(b(i,j+1,k)-lsv)/deltax;
        
        if(ymin*lsSig>0.0 && yplus*lsSig>-ymin*lsSig)
        dy=ddwenoy(a,b,1.0);

        if(yplus*lsSig<0.0 && ymin*lsSig<-yplus*lsSig)
        dy=ddwenoy(a,b,-1.0);

        if(yplus*lsSig>0.0 && ymin*lsSig<0.0)
        dy=0.0;

    // z
        zmin=(lsv-b(i,j,k-1))/deltax;
        zplus=(b(i,j,k+1)-lsv)/deltax;
        
        if(zmin*lsSig>0.0 && zplus*lsSig>-zmin*lsSig)
        dz=ddwenoz(a,b,1.0);

        if(zplus*lsSig<0.0 && zmin*lsSig<-zplus*lsSig)
        dz=ddwenoz(a,b,-1.0);

        if(zplus*lsSig>0.0 && zmin*lsSig<0.0)
        dz=0.0;

		
		dnorm = sqrt(dx*dx + dy*dy + dz*dz);


		lambda1 =  dirac * ((b(i,j,k) - d0(i,j,k))/dT) * dV;
		
		lambda2 =  dirac*dirac * dnorm * dV;
		
		
		lambda2 = fabs(lambda2)>1.0e-19?lambda2:1.0e10;
		
		b(i,j,k)      -=  p->F39*dT*dirac*dnorm*(lambda1/lambda2);
	}

	
}

void reinigc_RK3::step(fdm* a, lexer* p)
{
	dt = p->F43*deltax;
	dT = p->F43*deltax;
	reiniter=p->F44;
}

void reinigc_RK3::wallf_update(lexer *p, fdm *a, ghostcell *pgc)
{
	int n;
	FLUIDLOOP
	wallf(i,j,k)=0;
	
	GC4LOOP
	if(p->gcb4[n][4]==21 || p->gcb4[n][4]==22 || p->gcb4[n][4]==5)
	{
	i = p->gcb4[n][0];
	j = p->gcb4[n][1];
	k = p->gcb4[n][2];
	
	wallf(i,j,k)=1;
	}
}
