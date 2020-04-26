/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"levelset_RK3_V.h"
#include"gradient.h"
#include"lexer.h"
#include"field.h"
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
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS4.h"
#include"flux_HJ_CDS2_vrans.h"

levelset_RK3_V::levelset_RK3_V(lexer* p, fdm *a, ghostcell* pgc, heat *&pheat, concentration *&pconc):ddweno_nug(p),ark1(p),ark2(p),f(p),L(p)
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
    
    
    if(p->B269==0 && p->D11!=4)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269==0 && p->D11==4)
    pflux = new flux_HJ_CDS4(p);
    
    if(p->B269>=1 || p->S10==2)
    pflux = new flux_HJ_CDS2_vrans(p);
}

levelset_RK3_V::~levelset_RK3_V()
{
}

void levelset_RK3_V::start(fdm* a,lexer* p, convection* pconvec,solver* psolv, ghostcell* pgc,ioflow* pflow, reini* preini, particlecorr* ppart, field &ls)
{
    /*
    pflow->fsfinflow(p,a,pgc);
    pflow->fsfrkin(p,a,pgc,ark1);
    pflow->fsfrkin(p,a,pgc,ark2);
    pflow->fsfrkout(p,a,pgc,ark1);
    pflow->fsfrkout(p,a,pgc,ark2);
    ppicard->volcalc(p,a,pgc,ls);*/
	
	pflow->phi_relax(p,pgc,ls);

//  Fill V
    n=0;
	FLUIDLOOP
	{
	f.V[n]=ls(i,j,k);
	++n;
	}
    
	pgc->start4V(p,f,gcval_phi);


// Step 1
    starttime=pgc->timer();
    
    disc(p,a,f);

	NLOOP4
	ark1.V[n] = f.V[n]
				+ p->dt*L.V[n];
	
	//pflow->phi_relax(p,pgc,ark1);
	
    pgc->start4V(p,ark1,gcval_phi);
    
// Step 2
    disc(p,a,ark1);

	NLOOP4
	ark2.V[n] = 0.75*f.V[n]
				   + 0.25*ark1.V[n]
				   + 0.25*p->dt*L.V[n];
				
	//pflow->phi_relax(p,pgc,ark2);
	
	 pgc->start4V(p,ark2,gcval_phi);

// Step 3
    disc(p,a,ark2);

	NLOOP4
	f.V[n] =     (1.0/3.0)*f.V[n]
				  + (2.0/3.0)*ark2.V[n]
				  + (2.0/3.0)*p->dt*L.V[n];
                  
    if(p->count%p->F41==0)
	preini->startV(a,p,f,pgc,pflow);
                  
                  
// backfill
	n=0;
	FLUIDLOOP
	{
	ls(i,j,k)=f.V[n];
	++n;
	}

    pflow->phi_relax(p,pgc,ls);
	pgc->start4(p,ls,gcval_phi);

    ppart->start(p,a,pgc,pflow);
    pgc->start4(p,ls,gcval_phi);
	
	p->lsmtime=pgc->timer()-starttime;

    
	

    ppicard->correct_ls(p,a,pgc,ls);
	ppart->picardmove(p,a,pgc);

	pupdate->start(p,a,pgc);

	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"lsmtime: "<<setprecision(3)<<p->lsmtime<<endl;
}

void levelset_RK3_V::ltimesave(lexer* p, fdm *a, field &ls)
{
}

void levelset_RK3_V::update(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    pupdate->start(p,a,pgc);
}

void levelset_RK3_V::disc(lexer *p, fdm *a, vec &b)
{
	double dx,dy,dz;
    double iadvec,jadvec,kadvec;
    double ivel2,jvel2,kvel2;
    
n=0;
FLUIDLOOP
{
    pflux->u_flux(a,4,a->u,iadvec,ivel2);
    pflux->v_flux(a,4,a->v,jadvec,jvel2);
    pflux->w_flux(a,4,a->w,kadvec,kvel2);
        
	dx=0.0;
	dy=0.0;
	dz=0.0;

// x	
	if(iadvec>0.0)
	dx=iadvec*ddwenox(a,b,1.0,a->C4);

	if(iadvec<0.0)
	dx=iadvec*ddwenox(a,b,-1.0,a->C4);

// y
    if(p->j_dir==1)
    {
	if(jadvec>0.0)
	dy=jadvec*ddwenoy(a,b,1.0,a->C4);

	if(jadvec<0.0)
	dy=jadvec*ddwenoy(a,b,-1.0,a->C4);
    }

// z

	if(kadvec>0.0)
	dz=kadvec*ddwenoz(a,b,1.0,a->C4);

	if(kadvec<0.0)
	dz=kadvec*ddwenoz(a,b,-1.0,a->C4);


	L.V[n] = -dx-dy-dz;
    
    ++n;
}

}
