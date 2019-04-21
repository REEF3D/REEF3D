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

#include"kepsilon_AB.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"
#include"solver.h"
#include"diffusion.h"
#include"ioflow.h"
#include"kepsilon.h"
#include"convection.h"

kepsilon_AB::kepsilon_AB(lexer* p, fdm* a, ghostcell *pgc) : kepsilon(p,a,pgc),kab(p),eab(p)
{
	gcval_kin=20;
	gcval_eps=30;
	
	if(p->B67>0)
	gcval_kin=21;
	
	if(p->B68>0)
	gcval_eps=31;
}

kepsilon_AB::~kepsilon_AB()
{
}

void kepsilon_AB::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff,solver* psolv, ghostcell* pgc, ioflow* pflow)
{	
	Pk_update(p,a,pgc);
	wallf_update(p,a,pgc,wallf);

//kin
    starttime=pgc->timer();

	kinsource(p,a,kin,eps);
	bckeps_start(a,p,kin,eps,gcval_kin);
	pconvec->start(p,a,kin,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,kin,a->visc,ke_sigma_k,1.0);

	if(p->count==1)
	LOOP
	kab(i,j,k)=a->L(i,j,k);


	LOOP
	{
	kin(i,j,k)+=p->dt*0.5*(((p->turbtimestep+2.0*p->turbtimestep_old)/p->dt_old)*a->L(i,j,k)
								-(p->turbtimestep/p->turbtimestep_old)*kab(i,j,k));
	kab(i,j,k)=a->L(i,j,k);
	}
	
	pdiff->idiff_scalar(p,a,pgc,psolv,kin,a->visc,ke_sigma_k,1.0);
	pgc->start4(p,kin,gcval_kin);

	p->kintime=pgc->timer()-starttime;


//eps
    starttime=pgc->timer();

	epssource(p,a,kin,eps);
	pconvec->start(p,a,eps,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,eps,a->visc,ke_sigma_e,1.0);

	if(p->count==1)
	LOOP
	eab(i,j,k)=a->L(i,j,k);


	LOOP
	{
	eps(i,j,k)+=p->dt*0.5*(((p->turbtimestep+2.0*p->turbtimestep_old)/p->turbtimestep_old)*a->L(i,j,k)
								-(p->turbtimestep/p->turbtimestep_old)*eab(i,j,k));
	eab(i,j,k)=a->L(i,j,k);
	}
	
	pdiff->idiff_scalar(p,a,pgc,psolv,eps,a->visc,ke_sigma_e,1.0);
	epsfsf(p,a,pgc,kin,eps);
    bckeps_start(a,p,kin,eps,gcval_eps);
	pgc->start4(p,eps,gcval_eps);

	p->epstime=pgc->timer()-starttime;

	eddyvisc(a,p,kin,eps,pgc);
	
	pgc->start4(p,a->eddyv,24);
}

void kepsilon_AB::ktimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

void kepsilon_AB::etimesave(lexer *p, fdm* a, ghostcell *pgc)
{
}

