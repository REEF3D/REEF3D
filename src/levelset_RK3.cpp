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
#include"levelset_RK3.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reini.h"
#include"particlecorr.h"
#include"picard.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_void.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"

levelset_RK3::levelset_RK3(lexer* p, fdm *a, ghostcell* pgc, heat *&pheat, concentration *&pconc):gradient(p),ark1(p),ark2(p)
{
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;

	if(p->F30>0 && p->H10==0 && p->W30==0 && p->W90==0)
	pupdate = new fluid_update_fsf(p,a,pgc);
	
	if(p->F30>0 && p->H10==0 && p->W30==1 && p->W90==0)
	pupdate = new fluid_update_fsf_comp(p,a,pgc);
	
	if(p->F30>0 && p->H10>0 && p->W90==0 && p->H3==1)
	pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);
    
    if(p->F30>0 && p->H10>0 && p->W90==0 && p->H3==2)
	pupdate = new fluid_update_fsf_heat_Bouss(p,a,pgc,pheat);
	
	if(p->F30>0 && p->C10>0 && p->W90==0)
	pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);
	
	if(p->F30>0 && p->H10==0 && p->W30==0 && p->W90>0)
	pupdate = new fluid_update_rheology(p,a);

	if(p->F46==2)
	ppicard = new picard_f(p);

	if(p->F46==3)
	ppicard = new picard_lsm(p);

	if(p->F46!=2 && p->F46!=3)
	ppicard = new picard_void(p);
    
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    if(p->B22==1)
    {
    gcval_uls=117;
	gcval_vls=118;
	gcval_wls=119;
    }
    
    if(p->B22==2)
    {
    gcval_uls=110;
	gcval_vls=111;
	gcval_wls=112;
    }
}

levelset_RK3::~levelset_RK3()
{
}

void levelset_RK3::start(fdm* a,lexer* p, convection* pconvec,solver* psolv, ghostcell* pgc,ioflow* pflow, reini* preini, particlecorr* ppart, field &ls)
{
    pflow->fsfinflow(p,a,pgc);
    pflow->fsfrkin(p,a,pgc,ark1);
    pflow->fsfrkin(p,a,pgc,ark2);
    pflow->fsfrkout(p,a,pgc,ark1);
    pflow->fsfrkout(p,a,pgc,ark2);
    ppicard->volcalc(p,a,pgc,ls);
    
    pgc->start1(p,a->u,gcval_uls);
	pgc->start2(p,a->v,gcval_vls);
	pgc->start3(p,a->w,gcval_wls);
	

// Step 1
    starttime=pgc->timer();

    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pconvec->start(p,a,ls,4,a->u,a->v,a->w);
	

	FLUIDLOOP
	ark1(i,j,k) = ls(i,j,k)
				+ p->dt*a->L(i,j,k);
	
	pflow->phi_relax(p,pgc,ark1);
	
	pgc->start4(p,ark1,gcval_phi);
    
// Step 2
    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pconvec->start(p,a,ark1,4,a->u,a->v,a->w);

	FLUIDLOOP
	ark2(i,j,k) = 0.75*ls(i,j,k)
				   + 0.25*ark1(i,j,k)
				   + 0.25*p->dt*a->L(i,j,k);
				
	pflow->phi_relax(p,pgc,ark2);
	
	pgc->start4(p,ark2,gcval_phi);

// Step 3
    FLUIDLOOP
	a->L(i,j,k)=0.0;

	pconvec->start(p,a,ark2,4,a->u,a->v,a->w);

	FLUIDLOOP
	ls(i,j,k) =     (1.0/3.0)*ls(i,j,k)
				  + (2.0/3.0)*ark2(i,j,k)
				  + (2.0/3.0)*p->dt*a->L(i,j,k);

    pflow->phi_relax(p,pgc,ls);
	pgc->start4(p,ls,gcval_phi);

    ppart->start(p,a,pgc,pflow);
    
	
	p->lsmtime=pgc->timer()-starttime;

	preini->start(a,p,ls, pgc, pflow);
	

    ppicard->correct_ls(p,a,pgc,ls);
	ppart->picardmove(p,a,pgc);
    
	pupdate->start(p,a,pgc);

	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"lsmtime: "<<setprecision(3)<<p->lsmtime<<endl;
    
    pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);
}

void levelset_RK3::ltimesave(lexer* p, fdm *a, field &ls)
{
}

void levelset_RK3::update(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    pupdate->start(p,a,pgc);
}

