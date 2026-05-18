/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"nhflow_concentration_RK2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

nhflow_concentration_RK2::nhflow_concentration_RK2(lexer* p, fdm_nhf *d, ghostcell *pgc) : nhflow_concentration_io(p,d)
{
	gcval_concentration=80;
}

nhflow_concentration_RK2::~nhflow_concentration_RK2()
{
}

void nhflow_concentration_RK2::start(lexer *p, fdm_nhf *d, convection* pconvec, diffusion* pdiff, turbulence *pturb, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    /*field4 ark1(p);
    
    // Step 1
    starttime=pgc->timer();

    clearrhs(p,a,pgc);
    pconvec->start(p,a,C,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,C,a->visc,a->eddyv,1.0,1.0);

	LOOP
	ark1(i,j,k) = C(i,j,k)
				+ p->dt*a->L(i,j,k);

    bc_concentration_start(p,a,pgc,ark1);
	pgc->start4(p,ark1,gcval_concentration);


// Step 2
    clearrhs(p,a,pgc);
    pconvec->start(p,a,ark1,4,a->u,a->v,a->w);
	pdiff->diff_scalar(p,a,pgc,psolv,ark1,a->visc,a->eddyv,1.0,0.5);

	LOOP
	C(i,j,k) = 0.5*C(i,j,k)
				+ 0.5*ark1(i,j,k)
				+ 0.5*p->dt*a->L(i,j,k);
	
    bc_concentration_start(p,a,pgc,C);
	pgc->start4(p,C,gcval_concentration);

	pupdate->start(p,a,pgc,a->u,a->v,a->w);

	p->concentrationtime=pgc->timer()-starttime;*/

}

void nhflow_concentration_RK2::ctimesave(lexer *p, fdm_nhf *d)
{
}

void nhflow_concentration_RK2::clearrhs(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    /*
    int n=0;
	LOOP
	{
    a->L(i,j,k)=0.0;
	a->rhsvec.V[n]=0.0;
	++n;
	}*/
}
