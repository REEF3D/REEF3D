/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"levelset_RK2.h"
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
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"

levelset_RK2::levelset_RK2(lexer* p, fdm *a, ghostcell* pgc, heat *&pheat, concentration *&pconc):gradient(p)
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
}

levelset_RK2::~levelset_RK2()
{
}

void levelset_RK2::start(fdm* a,lexer* p, convection* pconvec,solver* psolv, ghostcell* pgc,ioflow* pflow, reini* preini, particlecorr* ppart, field &ls)
{
    field4 ark1(p);
    pflow->fsfinflow(p,a,pgc);
    pflow->fsfrkin(p,a,pgc,ark1);
    pflow->fsfrkout(p,a,pgc,ark1);
    ppicard->volcalc(p,a,pgc,ls);

// Step 1
    starttime=pgc->timer();

    LOOP
	a->L(i,j,k)=0.0;

	pconvec->start(p,a,ls,4,a->u,a->v,a->w);

	LOOP
	ark1(i,j,k) = ls(i,j,k)
				+ p->dt*a->L(i,j,k);

	pflow->phi_relax(p,pgc,ark1);
	
	pgc->start4(p,ark1,gcval_phi);

// Step 2
    LOOP
	a->L(i,j,k)=0.0;

	pconvec->start(p,a,ark1,4,a->u,a->v,a->w);

	LOOP
	ls(i,j,k) = 0.5*ls(i,j,k)
				  + 0.5*ark1(i,j,k)
				  + 0.5*p->dt*a->L(i,j,k);

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
}

void levelset_RK2::ltimesave(lexer* p, fdm *a, field &ls)
{
}

void levelset_RK2::update(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    pupdate->start(p,a,pgc);
}
