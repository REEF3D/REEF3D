/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"suspended_IM2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"sediment.h"

suspended_IM2::suspended_IM2(lexer* p, fdm* a, turbulence *pturb) : ibcsusp(p,pturb),isusprhs(p),concn(p),concnn(p),wvel(p)
{
	gcval_susp=60;
}

suspended_IM2::~suspended_IM2()
{
}

void suspended_IM2::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow, sediment *psed)
{
    starttime=pgc->timer();
    clearrhs(p,a);
    fill_wvel(p,a,pgc.psed);
    pconvec->start(p,a,a->conc,4,a->u,a->v,wvel);
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->eddyv,1.0,1.0);
	isuspsource(p,a,a->conc);
	timesource(p,a,a->conc);
	psolv->start(p,a,pgc,a->conc,a->rhsvec,4);
	ibcsusp_start(p,a,pgc,a->conc);
	sedfsf(p,a,a->conc);
	pgc->start4(p,a->conc,gcval_susp);
	p->susptime=pgc->timer()-starttime;
	p->suspiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0) && (p->count%p->P12==0))
	cout<<"suspiter: "<<p->suspiter<<"  susptime: "<<setprecision(3)<<p->susptime<<endl;
}

void suspended_IM2::timesource(lexer* p, fdm* a, field& fn)
{
    int count=0;
    int q;

    count=0;
    LOOP
    {
        a->M.p[count]+= 1.5/DT;

        a->rhsvec.V[count] += a->L(i,j,k) + (2.0*concn(i,j,k))/DT - concnn(i,j,k)/(2.0*DT);

	++count;
    }
}

void suspended_IM2::ctimesave(lexer *p, fdm* a)
{
    LOOP
    {
    concnn(i,j,k)=concn(i,j,k);
    concn(i,j,k)=a->conc(i,j,k);
    }

}

void suspended_IM2::fill_wvel(lexer *p, fdm* a, ghostcell *pgc, sediment *psed)
{
    WLOOP
    wvel(i,j,k) = a->w(i,j,k) + ws;
    
    pgc->start3(p,wvel,12);
}

