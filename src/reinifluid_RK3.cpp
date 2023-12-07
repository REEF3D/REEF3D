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

#include"reinifluid_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"
#include"reinidisc_f.h"
#include"reinidisc_sf.h"
#include"reinidisc_f2.h"
#include"reinidisc_fsf.h"

reinifluid_RK3::reinifluid_RK3(lexer* p, int type) : epsi(p->F45*p->DXM),f(p),frk1(p),frk2(p),L(p),dt(p)
{
	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=54;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

    gcval_iniphi=50;

    if((p->F61>1.0e-20 || p->F60>1.0e-20) && p->F50==1)
    gcval_iniphi=51;
	
	if((p->F62>1.0e-20 || p->F60>1.0e-20) && p->F50==2)
    gcval_iniphi=52;
	
	if(((p->F61>1.0e-20 && p->F62>1.0e-20) || p->F60>1.0e-20) && p->F50==3)
    gcval_iniphi=53;

	gcval_ro=1;

	if(type==41)
	gcval_iniphi=50;

    
	if(p->F46==1)
	ppicard = new picard_f(p);

	if(p->F46!=1)
	ppicard = new picard_void(p);
	
	if(p->F49==0)
	prdisc = new reinidisc_fsf(p);

	if(p->F49==1 && p->G3==0)
	prdisc = new reinidisc_f(p);
    
    if(p->F49==1 && p->G3==1)
	prdisc = new reinidisc_sf(p);
    
    if(p->F49==2)
	prdisc = new reinidisc_f2(p);
    
    
    time_preproc(p);
}

reinifluid_RK3::~reinifluid_RK3()
{
}

void reinifluid_RK3::start(fdm* a,lexer* p,field& b,ghostcell* pgc,ioflow* pflow)
{ 
	sizeM=p->sizeM4;
	
	ppicard->volcalc(p,a,pgc,a->phi);
	
	n=0;
	FLUIDLOOP
	{
	f.V[n]=b(i,j,k);
	++n;
	}
    
	pgc->start4vec(p,f,gcval_iniphi);
    
    startV(a,p,f,pgc,pflow);
    
    
    // backfill
	n=0;
	FLUIDLOOP
	{
	b(i,j,k)=f.V[n];
	++n;
	}
    
	if(p->count==0)
	pgc->start4(p,b,gcval_iniphi);
    
    if(p->count>0)
	pgc->start4(p,b,gcval_phi);
    
    ppicard->correct_ls(p,a,pgc,a->phi);
}

void reinifluid_RK3::startV(fdm* a,lexer* p,vec &f, ghostcell* pgc,ioflow* pflow)
{ 
    starttime=pgc->timer();
    
	if(p->count==0)
	{
    if(p->mpirank==0)
	cout<<"initializing level set..."<<endl<<endl;
	reiniter=2*int(p->maxlength/(p->F43*p->DXM));
	pgc->start4vec(p,f,gcval_iniphi);
	pflow->fsfrkinV(p,a,pgc,f);
	pflow->fsfrkoutV(p,a,pgc,f);
	}

	if(p->count>0)
	step(a,p);

	pflow->fsfrkinV(p,a,pgc,frk1);
    pflow->fsfrkinV(p,a,pgc,frk2);
    pflow->fsfrkoutV(p,a,pgc,frk1);
    pflow->fsfrkoutV(p,a,pgc,frk2);
    
    for(int q=0;q<reiniter;++q)
    {
	// Step 1
	prdisc->start(p,a,pgc,f,L,4);

	NLOOP4
	frk1.V[n] = f.V[n] + dt.V[n]*L.V[n];

	if(p->count==0)
	pgc->start4vec(p,frk1,gcval_iniphi);
    
    if(p->count>0)
	pgc->start4vec(p,frk1,gcval_phi);

    // Step 2
    prdisc->start(p,a,pgc,frk1,L,4);

	NLOOP4
	frk2.V[n] = 0.75*f.V[n] + 0.25*frk1.V[n] + 0.25*dt.V[n]*L.V[n];

	if(p->count==0)
	pgc->start4vec(p,frk2,gcval_iniphi);
    
    if(p->count>0)
	pgc->start4vec(p,frk2,gcval_phi);

    // Step 3
    prdisc->start(p,a,pgc,frk2,L,4);

	NLOOP4
	f.V[n] = (1.0/3.0)*f.V[n] + (2.0/3.0)*frk2.V[n] + (2.0/3.0)*dt.V[n]*L.V[n];

	if(p->count==0)
	pgc->start4vec(p,f,gcval_iniphi);
    
    if(p->count>0)
	pgc->start4vec(p,f,gcval_phi);
	}
    
    p->reinitime+=pgc->timer()-starttime;  
}

void reinifluid_RK3::step(fdm* a, lexer* p)
{
	reiniter=p->F44;
}

void reinifluid_RK3::time_preproc(lexer* p)
{
    n=0;
	FLUIDLOOP
	{
	dt.V[n] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	++n;
	}
}

