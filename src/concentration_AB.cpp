/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"concentration_AB.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_concentration.h"

concentration_AB::concentration_AB(lexer* p, fdm* a, ghostcell *pgc) : bc_concentration(p), concentration_io(p,a), cab(p)
{
	gcval_concentration=80;
}

concentration_AB::~concentration_AB()
{
}

void concentration_AB::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, turbulence *pturb, solver* psolv, ghostcell* pgc, ioflow* pflow)
{

    starttime=pgc->timer();
    
    clearrhs(p,a,pgc);
	pconvec->start(p,a,C,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,C,a->visc,a->eddyv,1.0,1.0);

	if(p->count==1)
	LOOP
	cab(i,j,k)=a->L(i,j,k);

	LOOP
	{
	C(i,j,k)+=p->dt*0.5*(((p->dt+2.0*p->dt_old)/p->dt_old)*a->L(i,j,k)
								-(p->dt/p->dt_old)*cab(i,j,k));
	cab(i,j,k)=a->L(i,j,k);
	}
	
    bc_concentration_start(p,a,pgc,C);
	pgc->start4(p,C,gcval_concentration);

	pupdate->start(p,a,pgc);

	p->concentrationtime=pgc->timer()-starttime;
}

void concentration_AB::ttimesave(lexer *p, fdm* a)
{
}

void concentration_AB::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}
}
