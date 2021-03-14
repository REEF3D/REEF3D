/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"suspended_RK3.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

suspended_RK3::suspended_RK3(lexer* p, fdm* a, turbulence *pturb) : bcsusp(p,pturb),susprhs(p)
{
	gcval_susp=60;
}

suspended_RK3::~suspended_RK3()
{
}


void suspended_RK3::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    field4 ark1(p),ark2(p);
    
// Step 1
    starttime=pgc->timer();

    suspsource(p,a,a->conc);
    pconvec->start(p,a,a->conc,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,a->conc,a->visc,1.0,1.0);

	LOOP
	ark1(i,j,k) = a->conc(i,j,k)
                   + p->dt*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,ark1,a->visc,1.0,1.0);
    bcsusp_start(p,a,pgc,ark1);
    sedfsf(p,a,ark1);
	pgc->start4(p,ark1,gcval_susp);

// Step 2
    suspsource(p,a,a->conc);
    pconvec->start(p,a,ark1,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,ark1,a->visc,1.0,0.25);

	LOOP
	ark2(i,j,k) = 0.75*a->conc(i,j,k)
                   + 0.25*ark1(i,j,k)
				   + 0.25*p->dt*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,ark2,a->visc,1.0,0.25);
    bcsusp_start(p,a,pgc,ark2);
    sedfsf(p,a,ark2);
	pgc->start4(p,ark2,gcval_susp);

// Step 3
    suspsource(p,a,a->conc);
    pconvec->start(p,a,ark2,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,ark2,a->visc,1.0,2.0/3.0);

	LOOP
	a->conc(i,j,k) = (1.0/3.0)*a->conc(i,j,k)
				+ (2.0/3.0)*ark2(i,j,k)
				+ (2.0/3.0)*p->dt*a->L(i,j,k);
	
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->visc,1.0,2.0/3.0);
    bcsusp_start(p,a,pgc,a->conc);
    sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);

	p->susptime=pgc->timer()-starttime;

}

void suspended_RK3::ctimesave(lexer *p, fdm* a)
{
}
