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

#include"suspended_IM2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"discrete.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

suspended_IM2::suspended_IM2(lexer* p, fdm* a, turbulence *pturb) : ibcsusp(p,pturb),isusprhs(p),concn(p),concnn(p)
{
	gcval_susp=60;
}

suspended_IM2::~suspended_IM2()
{
}

void suspended_IM2::start(fdm* a, lexer* p, discrete* pdisc, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
    starttime=pgc->timer();
    clearrhs(p,a);
    pdisc->start(p,a,a->conc,4,a->u,a->v,a->w);
	pdiff->idiff_scalar(p,a,pgc,psolv,a->conc,a->visc,1.0,1.0);
	isuspsource(p,a,a->conc);
	timesource(p,a,a->conc);
	psolv->start(p,a,pgc,a->conc,a->xvec,a->rhsvec,4,gcval_susp,p->N43);
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
        a->M.p[count]+= 1.5/PDT;

        a->rhsvec.V[count] += a->L(i,j,k) + (2.0*concn(i,j,k))/PDT - concnn(i,j,k)/(2.0*PDT);

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


