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

#include"reini_AB3.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"

reini_AB3::reini_AB3(lexer* p, fdm *a):f(p),dab1(p),dab2(p),L(p),dt(p)
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
	gcval_ro=1;

	
	if((p->F61>1.0e-20 || p->F60>1.0e-20) && p->F50==1)
    gcval_iniphi=51;
	
	if((p->F62>1.0e-20 || p->F60>1.0e-20) && p->F50==2)
    gcval_iniphi=52;

	if(p->F46==1)
	ppicard = new picard_f(p);

	if(p->F46!=1)
	ppicard = new picard_void(p);
	
	if(p->F49==0)
	prdisc = new reinidisc_fsf(p);

	if(p->F49==1)
	prdisc = new reinidisc_f(p);
    
    time_preproc(p);
}

reini_AB3::~reini_AB3()
{
}

void reini_AB3::start(fdm* a,lexer* p,field& b, ghostcell* pgc,ioflow* pflow)
{
	starttime=pgc->timer();
	
	sizeM=p->sizeM4;
	
	ppicard->volcalc(p,a,pgc,a->phi);
	
	 // fill lsm to reini
	n=0;
	BASELOOP
	{
	PFLUIDCHECK
	f.V[n]=b(i,j,k);
	
	++n;
	}

	pgc->start6V(p,f,gcval_iniphi);
	
	if(p->count==0)
	{
	inisolid(p,a);
    if(p->mpirank==0)
	cout<<"initializing level set..."<<endl<<endl;
	reiniter=2*int(p->maxlength/(p->F43*p->DXM));
	pgc->start6V(p,f,gcval_iniphi);
	fsfrkioV(p,a,pgc,f);
	}

	if(p->count>0)
	step(p,a);

	ppicard->volcalc(p,a,pgc,a->phi);

	for(int q=0;q<reiniter;++q)
	{

		prdisc->start(p,a,pgc,f,L,4);

		if(q==0)
		NLOOP6
		{
		dab1.V[n]=L.V[n];
		dab2.V[n]=L.V[n];
		}
		
		NLOOP6
		{		
		f.V[n] +=dt.V[n]*(1.0/12.0)*(23.0*L.V[n] - 16.0*dab1.V[n] + 5.0*dab2.V[n]);
		
		dab2.V[n]=dab1.V[n];
		dab1.V[n]=L.V[n];
		
		}

	if(p->count==0)
	pgc->start6V(p,f,gcval_iniphi);

    if(p->count>0)
	pgc->start6V(p,f,gcval_phi);
	}

    // backfill
	n=0;
	LOOP
	{
	b(i,j,k)=f.V[n];
	++n;
	}
	
	if(p->count==0)
	pgc->start4(p,b,gcval_iniphi);
    
    if(p->count>0)
	pgc->start4(p,b,gcval_phi);
	
	ppicard->correct_ls(p,a,pgc,a->phi);
	
	p->reinitime+=pgc->timer()-starttime;

}

void reini_AB3::step(lexer *p, fdm *a)
{
	reiniter=p->F44+2;
}

void reini_AB3::time_preproc(lexer* p)
{
	/*
	n=0;
	LOOP
	{
	dt.V[n] = p->F43*(1.0/3.0)*(p->DXP[IP] + p->DYP[JP] + p->DZP[KP]);
	++n;
	}*/
    
    n=0;
	LOOP
	{
	dt.V[n] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	++n;
	}
}


void reini_AB3::inisolid(lexer* p, fdm *a)
{
	pip=4;
	n=0;
	BASELOOP
	{
    f.V[n]=a->phi(i,j,k);
	 
	++n;
	}
	pip=0;
}

void reini_AB3::fsfrkioV(lexer *p, fdm *a, ghostcell *pgc, vec& f)
{
    
    for(int q=0;q<p->gcin6_count;++q)
    {
        i=p->gcin6[q][0];
        j=p->gcin6[q][1];
        k=p->gcin6[q][2];
        n=p->gcin6[q][3];

        f.V[Im1_J_K_6]=a->phi(i-1,j,k);
        f.V[Im2_J_K_6]=a->phi(i-2,j,k);
        f.V[Im3_J_K_6]=a->phi(i-3,j,k);
    }
    
    
    for(int q=0;q<p->gcout6_count;++q)
    {
        i=p->gcout6[q][0];
        j=p->gcout6[q][1];
        k=p->gcout6[q][2];
        n=p->gcout6[q][3];

        f.V[Ip1_J_K_6]=a->phi(i+1,j,k);
        f.V[Ip2_J_K_6]=a->phi(i+2,j,k);
        f.V[Ip3_J_K_6]=a->phi(i+3,j,k);
    }
}



