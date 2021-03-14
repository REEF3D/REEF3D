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

#include"suspended_AB.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

suspended_AB::suspended_AB(lexer* p, fdm* a, turbulence *pturb) : bcsusp(p,pturb),susprhs(p),cab(p)
{
	gcval_susp=60;
}

suspended_AB::~suspended_AB()
{
}

void suspended_AB::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{

    starttime=pgc->timer();

	suspsource(p,a,a->conc);
	pconvec->start(p,a,a->conc,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,a->conc,a->visc,1.0,1.0);

	if(p->count==1)
	LOOP
	cab(i,j,k)=a->L(i,j,k);


	LOOP
	{
	a->conc(i,j,k)+=p->dt*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->L(i,j,k)
								-(p->dt/p->dt_old)*cab(i,j,k));
	cab(i,j,k)=a->L(i,j,k);
	}
	
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->visc,1.0,1.0);
    bcsusp_start(p,a,pgc,a->conc);
    sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);

	p->susptime=pgc->timer()-starttime;
}

void suspended_AB::ctimesave(lexer *p, fdm* a)
{
}

