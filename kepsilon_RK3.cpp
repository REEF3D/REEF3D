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

#include"kepsilon_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"solver.h"
#include"diffusion.h"
#include"ioflow.h"
#include"kepsilon.h"
#include"convection.h"

kepsilon_RK3::kepsilon_RK3(lexer* p, fdm* a, ghostcell *pgc) : kepsilon(p,a,pgc)
{
	gcval_kin=20;
	gcval_eps=30;
	
	if(p->B67>0)
	gcval_kin=21;
	
	if(p->B68>0)
	gcval_eps=31;

	gcval_ark=21;
	gcval_brk=31;
}

kepsilon_RK3::~kepsilon_RK3()
{
}

void kepsilon_RK3::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow)
{	
    field4 ark1(p),ark2(p);
    field4 brk1(p),brk2(p);
	Pk_update(p,a,pgc);
	wallf_update(p,a,pgc,wallf);

// Step 1
	//kin
	starttime=pgc->timer();

	kinsource(p,a,kin,eps);
	bckeps_start(a,p,kin,eps,gcval_kin);
	pconvec->start(p,a,kin,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,kin,a->visc,ke_sigma_k,1.0);

	LOOP
	ark1(i,j,k) = kin(i,j,k)
				+ p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,ark1,a->visc,ke_sigma_k,1.0);
	pgc->start4(p,ark1,gcval_ark);

	p->kintime=pgc->timer()-starttime;

	//eps
	starttime=pgc->timer();

	epssource(p,a,kin,eps);
	pconvec->start(p,a,eps,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,eps,a->visc,ke_sigma_e,1.0);

	LOOP
	brk1(i,j,k) = eps(i,j,k)
				+ p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,brk1,a->visc,ke_sigma_e,1.0);
	epsfsf(p,a,pgc,brk1,ark1);
	bckeps_start(a,p,kin,eps,gcval_eps);
	pgc->start4(p,brk1,gcval_brk);

	p->epstime=pgc->timer()-starttime;

// Step 2
	//kin
	starttime=pgc->timer();

	kinsource(p,a,kin,eps);
	bckeps_start(a,p,kin,eps,gcval_kin);
	pconvec->start(p,a,ark1,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,ark1,a->visc,ke_sigma_k,0.25);

	LOOP
	ark2(i,j,k) = 0.75*kin(i,j,k) + 0.25*ark1(i,j,k)
				+ 0.25*p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,ark2,a->visc,ke_sigma_k,0.25);
	pgc->start4(p,ark2,gcval_ark);

	p->kintime+=pgc->timer()-starttime;

	//eps
	starttime=pgc->timer();

	epssource(p,a,kin,eps);
	pconvec->start(p,a,brk1,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,brk1,a->visc,ke_sigma_e,0.25);

	LOOP
	brk2(i,j,k) = 0.75*eps(i,j,k) + 0.25*brk1(i,j,k)
				+ 0.25*p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,brk2,a->visc,ke_sigma_e,0.25);
	epsfsf(p,a,pgc,brk2,ark2);
	bckeps_start(a,p,kin,eps,gcval_eps);
	pgc->start4(p,brk2,gcval_brk);

	p->epstime+=pgc->timer()-starttime;

// Step 3
	//kin
	starttime=pgc->timer();

	kinsource(p,a,kin,eps);
	bckeps_start(a,p,kin,eps,gcval_kin);
	pconvec->start(p,a,ark2,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,ark2,a->visc,ke_sigma_k,2.0/3.0);

	LOOP
	kin(i,j,k) = (1.0/3.0)*kin(i,j,k) + (2.0/3.0)*ark2(i,j,k)
				  + (2.0/3.0)*p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,kin,a->visc,ke_sigma_k,2.0/3.0);
	pgc->start4(p,kin,gcval_kin);

	p->kintime+=pgc->timer()-starttime;

	//eps
	starttime=pgc->timer();

	epssource(p,a,kin,eps);
	pconvec->start(p,a,brk2,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,brk2,a->visc,ke_sigma_e,2.0/3.0);

	LOOP
	eps(i,j,k) = (1.0/3.0)*eps(i,j,k) + (2.0/3.0)*brk2(i,j,k)
				  + (2.0/3.0)*p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,eps,a->visc,ke_sigma_e,2.0/3.0);
	epsfsf(p,a,pgc,eps,kin);
	bckeps_start(a,p,kin,eps,gcval_eps);
	pgc->start4(p,eps,gcval_eps);

	p->epstime+=pgc->timer()-starttime;

	eddyvisc(a,p,kin,eps,pgc);

	pflow->periodic(kin,p);
	pflow->periodic(eps,p);
	pgc->start4(p,a->eddyv,24);
}

void kepsilon_RK3::ktimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

void kepsilon_RK3::etimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

