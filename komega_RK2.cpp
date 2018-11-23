/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"komega_RK2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"solver.h"
#include"diffusion.h"
#include"ioflow.h"
#include"komega.h"
#include"convection.h"

komega_RK2::komega_RK2(lexer* p, fdm* a, ghostcell *pgc) : komega(p,a,pgc)
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

komega_RK2::~komega_RK2()
{
}

void komega_RK2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    field4 ark1(p);
    field4 brk1(p);
	Pk_update(p,a,pgc);
	wallf_update(p,a,pgc,wallf);

// Step 1
	//kin
	starttime=pgc->timer();

	kinsource(p,a,pgc,kin,eps);
    bckeps_start(a,p,kin,eps,gcval_kin);
    pconvec->start(p,a,kin,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,kin,a->visc,kw_sigma_k,1.0);

	LOOP
	ark1(i,j,k) = kin(i,j,k)
				+ p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,ark1,a->visc,kw_sigma_k,1.0);
	pgc->start4(p,ark1,gcval_ark);

	p->kintime=pgc->timer()-starttime;

	//omega
	starttime=pgc->timer();

	epssource(p,a,pgc,kin,eps);
	pconvec->start(p,a,eps,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,eps,a->visc,kw_sigma_w,1.0);

	LOOP
	brk1(i,j,k) = eps(i,j,k)
				+ p->turbtimestep*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,brk1,a->visc,kw_sigma_w,1.0);
	epsfsf(p,a,pgc,kin,eps);
	bckeps_start(a,p,kin,eps,gcval_eps);
	pgc->start4(p,brk1,gcval_brk);

	p->epstime=pgc->timer()-starttime;

// Step 2
	//kin
	starttime=pgc->timer();

    kinsource(p,a,pgc,kin,eps);
    bckeps_start(a,p,kin,eps,gcval_kin);
    pconvec->start(p,a,ark1,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,ark1,a->visc,kw_sigma_k,0.5);

	LOOP
	kin(i,j,k) = 0.5*kin(i,j,k) + 0.5*ark1(i,j,k)
				+ 0.5*p->turbtimestep*a->L(i,j,k);

	pdiff->idiff_scalar(p,a,pgc,psolv,kin,a->visc,kw_sigma_k,0.5);
	pgc->start4(p,kin,gcval_kin);

	p->kintime+=pgc->timer()-starttime;

	//omega
	starttime=pgc->timer();

	epssource(p,a,pgc,kin,eps);
	pconvec->start(p,a,brk1,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,brk1,a->visc,kw_sigma_w,0.5);

	LOOP
	eps(i,j,k) = 0.5*eps(i,j,k) + 0.5*brk1(i,j,k)
				+ 0.5*p->turbtimestep*a->L(i,j,k);

	pdiff->idiff_scalar(p,a,pgc,psolv,eps,a->visc,kw_sigma_w,0.5);
	epsfsf(p,a,pgc,kin,eps);
	bckeps_start(a,p,kin,eps,gcval_eps);
	pgc->start4(p,eps,gcval_eps);

	p->epstime+=pgc->timer()-starttime;

	eddyvisc(p,a,kin,eps,pgc);

	pflow->periodic(kin,p);
	pflow->periodic(eps,p);
	pgc->start4(p,a->eddyv,24);
}

void komega_RK2::ktimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

void komega_RK2::etimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

