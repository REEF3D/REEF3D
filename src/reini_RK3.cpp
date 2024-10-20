/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public a->License as published by
the Free Software Foundation; either version 3 of the a->License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABIa->LITY or
FITNESS FOR A PARTICUa->LAR PURPOSE. See the GNU General Public a->License
for more details.

You should have received a copy of the GNU General Public a->License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"reini_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"picard_f.h"
#include"picard_void.h"
#include"reinidisc_f.h"
#include"reinidisc_fsf.h"

reini_RK3::reini_RK3(lexer* p, int type) : epsi(p->F45*p->DXM),frk1(p),frk2(p),dt(p)
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

	if(p->F49==1)
	prdisc = new reinidisc_f(p);
    
    time_preproc(p);    
}

reini_RK3::~reini_RK3()
{
}

void reini_RK3::start(fdm *a, lexer *p, field &f, ghostcell *pgc, ioflow* pflow)
{ 
    starttime=pgc->timer();

    if(p->count==0)
    gcval = gcval_iniphi;
        
    if(p->count>0)
    gcval = gcval_phi;
    
	ppicard->volcalc(p,a,pgc,a->phi);
	
	if(p->count==0)
	{
        if(p->mpirank==0)
        cout<<"initializing level set..."<<endl<<endl;
        reiniter=2*int(p->maxlength/(p->F43*p->DXM));
        pgc->start4(p,f,gcval_iniphi);
	}
    
	if(p->count>0)
	step(p,a);
	
	pflow->fsfrkin(p,a,pgc,frk1);
    pflow->fsfrkin(p,a,pgc,frk2);
    pflow->fsfrkout(p,a,pgc,frk1);
    pflow->fsfrkout(p,a,pgc,frk2);

    
    for(int q=0;q<reiniter;++q)
    {
        // Step 1
        prdisc->start(p,a,pgc,f,a->L,4);

        BASELOOP
        frk1.V[IJK] = f.V[IJK] + dt.V[IJK]*a->L.V[IJK];
        
        pgc->start4(p,frk1,gcval);

        // Step 2
        prdisc->start(p,a,pgc,frk1,a->L,4);

        BASELOOP
        frk2.V[IJK] = 0.75*f.V[IJK] + 0.25*frk1.V[IJK] + 0.25*dt.V[IJK]*a->L.V[IJK];

        pgc->start4(p,frk2,gcval);

        // Step 3
        prdisc->start(p,a,pgc,frk2,a->L,4);

        BASELOOP
        f.V[IJK] = (1.0/3.0)*f.V[IJK] + (2.0/3.0)*frk2.V[IJK] + (2.0/3.0)*dt.V[IJK]*a->L.V[IJK];

        pgc->start4(p,f,gcval);
	}
    
	ppicard->correct_ls(p,a,pgc,a->phi);
	
	p->reinitime+=pgc->timer()-starttime;
}

void reini_RK3::step(lexer* p, fdm *a)
{
	reiniter=p->F44;
}

void reini_RK3::time_preproc(lexer* p)
{
	LOOP
	{
    if(p->j_dir==0)
    dt.V[IJK] = p->F43*MIN(p->DXP[IP],p->DZP[KP]);
    
    if(p->j_dir==1)
	dt.V[IJK] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	}
}


