/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

#include"nhflow_pjm_phs.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson_phs.h"
#include"density_f.h"
#include"patchBC_interface.h"

nhflow_pjm_phs::nhflow_pjm_phs(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(1.0)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson_phs(p);

    gcval_press=540;
    
    gamma=0.5;
    
    
    if(p->D33==0)
    solver_id = 8;
    
    if(p->D33==1)
    solver_id = 9;
    
}

nhflow_pjm_phs::~nhflow_pjm_phs()
{
}

void nhflow_pjm_phs::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL,
                        double *UH, double *VH, double *WH, double alpha)
{
    /*
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    rhs(p,d,pgc,d->U,d->V,d->W,alpha);

    ppois->start(p,d,d->PHS);

        starttime=pgc->timer();

    psolv->startF(p,pgc,d->PHS,d->rhsvec,d->M,solver_id);

        endtime=pgc->timer();

	pgc->start7P(p,d->PHS,gcval_press);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"phsiter: "<<p->solveriter<<"  phstime: "<<setprecision(3)<<p->poissontime<<endl;
    */
    
    FLOOP
    d->PHS[FIJK] = (p->wd + d->eta(i,j) - p->ZSN[FIJK])*p->W1*fabs(p->W22);
    
    pgc->start7P(p,d->PHS,gcval_press);
    
}

void nhflow_pjm_phs::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *UH, double *P, double alpha)
{
}

void nhflow_pjm_phs::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *VH, double *P, double alpha)
{
}

void nhflow_pjm_phs::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *WH, double *P, double alpha)
{
}

void nhflow_pjm_phs::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    n=0;
    FLOOP
    {
    d->rhsvec.V[n] =      -  p->W22;                  
    ++n;
    }
    
    //if(p->count<5)
    FLOOP
    d->PHS[FIJK] = (p->wd + d->eta(i,j) - p->ZSN[FIJK])*p->W1*fabs(p->W22);
    
    pgc->start7P(p,d->PHS,gcval_press);
}

void nhflow_pjm_phs::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_phs::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_phs::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
}

void nhflow_pjm_phs::vpgrad(lexer*p,fdm_nhf *d, slice &WL)
{
}

void nhflow_pjm_phs::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
}

