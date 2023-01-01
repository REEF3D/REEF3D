/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_void.h"
#include"heat.h"
#include"concentration.h"
#include"momentum.h"
#include"sflow_hxy_weno.h"
#include"sflow_hxy_cds.h"
#include"sflow_hxy_fou.h"
#include"patchBC_interface.h"

nhflow_fsf_rk::nhflow_fsf_rk(lexer *p, fdm_nhf* d, ghostcell *pgc, ioflow *pflow, patchBC_interface *ppBC) : epsi(p->A440*p->DXM),P(p),Q(p)
{
    pBC = ppBC;
    
	pupdate = new fluid_update_void();
    
    if(p->A541==1)
	phxy = new sflow_hxy_fou(p,pBC);
	
	if(p->A541==2)
	phxy = new sflow_hxy_cds(p,pBC);
	
	if(p->A541==4)
	phxy = new sflow_hxy_weno(p,pBC);
}

nhflow_fsf_rk::~nhflow_fsf_rk()
{
}

void nhflow_fsf_rk::start(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow)
{
}

void nhflow_fsf_rk::update(lexer *p, fdm_nhf* d, ghostcell *pgc, slice &f)
{
    //pupdate->start(p,a,pgc);
}

void nhflow_fsf_rk::ini(lexer *p, fdm_nhf* d, ghostcell *pgc, ioflow *pflow)
{
}

void nhflow_fsf_rk::step1(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += U[IJK]*p->DZN[KP];

    VLOOP
	Q(i,j) += V[IJK]*p->DZN[KP];
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,d->eta,P,Q);
    
    SLICELOOP1
    P(i,j) *= d->hx(i,j);

    SLICELOOP2
	Q(i,j) *= d->hy(i,j);

	pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);

    SLICELOOP4
    etark1(i,j) = d->eta(i,j)

                - p->dt*((P(i,j)-P(i-1,j))/p->DXN[IP] + (Q(i,j)-Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,etark1);
    pgc->gcsl_start4(p,etark1,1);
    
    //p->sigma_update(p,d,pgc,etark1,d->eta,1.0);
    //p->omega_update(p,d,pgc,u,v,w,etark1,d->eta,1.0);
}

void nhflow_fsf_rk::step2(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += U[IJK]*p->DZN[KP];

    VLOOP
	Q(i,j) += V[IJK]*p->DZN[KP];
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,etark1,P,Q);
    
    SLICELOOP1
    P(i,j) *= d->hx(i,j);

    SLICELOOP2
	Q(i,j) *= d->hy(i,j);

	pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);

    SLICELOOP4
    etark2(i,j) = 0.75*d->eta(i,j) + 0.25*etark1(i,j)

                - 0.25*p->dt*((P(i,j)-P(i-1,j))/p->DXN[IP] + (Q(i,j)-Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,etark2);
    pgc->gcsl_start4(p,etark2,1);
    
    //p->sigma_update(p,a,pgc,etark2,etark1,0.25);
    //p->omega_update(p,a,pgc,u,v,w,etark2,etark1,0.25);
}

void nhflow_fsf_rk::step3(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W, slice& etark1, slice &etark2, double alpha)
{
    // fill eta_n
    SLICELOOP4
    {
    d->eta_n(i,j) = d->eta(i,j);
    d->WL_n(i,j) = d->WL(i,j);
    }
    
    pgc->gcsl_start4(p,d->eta_n,1);
    
    SLICELOOP1
    P(i,j)=0.0;
    
    SLICELOOP2
    Q(i,j)=0.0;

    ULOOP
    P(i,j) += U[IJK]*p->DZN[KP];

    VLOOP
	Q(i,j) += V[IJK]*p->DZN[KP];
    
    pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);
    
    phxy->start(p,d->hx,d->hy,d->depth,p->wet,etark2,P,Q);
    
    SLICELOOP1
    P(i,j) *= d->hx(i,j);

    SLICELOOP2
	Q(i,j) *= d->hy(i,j);

	pgc->gcsl_start1(p,P,10);
    pgc->gcsl_start2(p,Q,11);

    SLICELOOP4
    d->eta(i,j) = (1.0-alpha)*d->eta(i,j) + alpha*etark2(i,j)

                - alpha*p->dt*((P(i,j)-P(i-1,j))/p->DXN[IP] + (Q(i,j)-Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,d->eta);
    pgc->gcsl_start4(p,d->eta,1);
    
    SLICELOOP4
    d->WL(i,j) = MAX(0.0, d->eta(i,j) + p->wd - d->bed(i,j));
    
    //p->sigma_update(p,a,pgc,d->eta,etark2,alpha);
    //p->omega_update(p,a,pgc,u,v,w,d->eta,etark2,alpha);
}


