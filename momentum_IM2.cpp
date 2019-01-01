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

#include"momentum_IM2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

momentum_IM2::momentum_IM2(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow)
                                                    :ibcmom(p),un(p),unn(p),vn(p),vnn(p),wn(p),wnn(p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

    utimesave(p,a,pgc);
	vtimesave(p,a,pgc);
	wtimesave(p,a,pgc);

	pconvec=pconvection;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
}

momentum_IM2::~momentum_IM2()
{
}

void momentum_IM2::start(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom)
{
    pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	
	p->uiter=p->viter=p->witer=0;
	p->utime=p->vtime=p->wtime=0.0;

	//--------------------------------------------------------
	//U
	starttime=pgc->timer();
	clearrhs(p,a);
	pturb->isource(p,a);
	pconvec->start(p,a,a->u,1,a->u,a->v,a->w);
	pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
	pflow->isource(p,a,pgc);
    ibcmom_start(a,p,pgc,pturb,a->u,gcval_u);
    ppress->upgrad(p,a);
    usource(p,a);
    ppress->fillapu(p,a);
	psolv->start(p,a,pgc,a->u,a->xvec,a->rhsvec,1,gcval_u,p->N43);
	pgc->start1(p,a->u,gcval_u);
	p->utime+=pgc->timer()-starttime;
	p->uiter+=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"uiter: "<<p->uiter<<"  utime: "<<setprecision(3)<<p->utime<<endl;

	//--------------------------------------------------------
	//V
	starttime=pgc->timer();
	clearrhs(p,a);
	pturb->jsource(p,a);
	pconvec->start(p,a,a->v,2,a->u,a->v,a->w);
	pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
	pflow->jsource(p,a,pgc);
    ibcmom_start(a,p,pgc,pturb,a->v,gcval_v);
    ppress->vpgrad(p,a);
    vsource(p,a);
    ppress->fillapv(p,a);
	psolv->start(p,a,pgc,a->v,a->xvec,a->rhsvec,2,gcval_v,p->N43);
	pgc->start2(p,a->v,gcval_v);
	p->vtime+=pgc->timer()-starttime;
	p->viter+=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"viter: "<<p->viter<<"  vtime: "<<setprecision(3)<<p->vtime<<endl;

	//--------------------------------------------------------
	//W
	starttime=pgc->timer();
	clearrhs(p,a);
	pturb->ksource(p,a);
	pconvec->start(p,a,a->w,3,a->u,a->v,a->w);
	pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,1.0);
	pflow->ksource(p,a,pgc);
    ibcmom_start(a,p,pgc,pturb,a->w,gcval_w);
    ppress->wpgrad(p,a);
    wsource(p,a);
    ppress->fillapw(p,a);
	psolv->start(p,a,pgc,a->w,a->xvec,a->rhsvec,3,gcval_w,p->N43);
	pgc->start3(p,a->w,gcval_w);
	p->wtime+=pgc->timer()-starttime;
	p->witer+=p->solveriter;
	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"witer: "<<p->witer<<"  wtime: "<<setprecision(3)<<p->wtime<<endl;

	//--------------------------------------------------------
	// pressure
	pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,1.0);

	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	pflow->periodic(a->u,p);
	pflow->periodic(a->v,p);
	pflow->periodic(a->w,p);

}




