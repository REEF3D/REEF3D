/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"concentration_IM1.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_fsf_concentration.h"

concentration_IM1::concentration_IM1(lexer* p, fdm* a, ghostcell *pgc) : ibc_concentration(p), concentration_io(p,a), Cn(p)
{
	gcval_concentration=80;
}

concentration_IM1::~concentration_IM1()
{
}

void concentration_IM1::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, turbulence *pturb, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    clearrhs(p,a,pgc);
    pconvec->start(p,a,C,4,a->u,a->v,a->w);
	pdiff->idiff_scalar(p,a,pgc,psolv,C,a->visc,1.0,1.0);
	timesource(p,a,C);
	psolv->start(p,a,pgc,C,a->xvec,a->rhsvec,4,gcval_concentration,p->N43);
	ibc_concentration_start(p,a,pgc,C);
	pgc->start4(p,C,gcval_concentration);
	pupdate->start(p,a,pgc);
	p->concentrationtime=pgc->timer()-starttime;
	p->concentrationiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"concentrationiter: "<<p->suspiter<<"  concentrationtime: "<<setprecision(3)<<p->susptime<<endl;
}

void concentration_IM1::timesource(lexer* p, fdm* a, field& fn)
{
    int count=0;

    LOOP
    {
        a->M.p[count]+= 1.0/PDT;

        a->rhsvec.V[count] += a->L(i,j,k) + Cn(i,j,k)/PDT;

	++count;
    }
}

void concentration_IM1::ttimesave(lexer* p, fdm* a)
{
    LOOP
    Cn(i,j,k)=C(i,j,k);
}

void concentration_IM1::clearrhs(lexer *p, fdm *a, ghostcell *pgc)
{

    int count=0;
    LOOP
    {
    a->L(i,j,k)=0.0;
    a->rhsvec.V[count]=0.0;
	++count;
    }

}

