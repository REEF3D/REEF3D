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

#include"nhflow_fsf_fsm.h"
#include"lexer.h"
#include"fdm.h"
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

nhflow_fsf_fsm::nhflow_fsf_fsm(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, patchBC_interface *ppBC) : epsi(p->A440*p->DXM), hx(p), hy(p)
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

nhflow_fsf_fsm::~nhflow_fsf_fsm()
{
}

void nhflow_fsf_fsm::start(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow)
{
    pgc->start4(p,a->ro,1);
    pgc->start4(p,a->visc,1);
    
    // fill eta_n
    SLICELOOP4
    a->eta_n(i,j) = a->eta(i,j);

    pgc->gcsl_start4(p,a->eta_n,gcval_phi);
    
    
    // Calculate Eta
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += a->u(i,j,k)*p->DZN[KP];

    VLOOP
	a->Q(i,j) += a->v(i,j,k)*p->DZN[KP];
    
    pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
    
    phxy->start(p,hx,hy,a->depth,a->wet,a->eta,a->P,a->Q);
    
    SLICELOOP1
    a->P(i,j) *= hx(i,j);

    SLICELOOP2
	a->Q(i,j) *= hy(i,j);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    // fsf equation
    SLICELOOP4
    a->eta(i,j) -= p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);	  
    
    pflow->eta_relax(p,pgc,a->eta);
    pgc->gcsl_start4(p,a->eta,1);
    
    p->sigma_update(p,a,pgc,a->eta,1.0);
    p->omega_update(p,a,pgc,a->u,a->v,a->w);
    
}

void nhflow_fsf_fsm::ltimesave(lexer* p, fdm *a, slice &ls)
{
}

void nhflow_fsf_fsm::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
    pupdate->start(p,a,pgc);
}

void nhflow_fsf_fsm::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}

void nhflow_fsf_fsm::step1(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, slice& etark1, slice &etark2, double alpha)
{
    SLICELOOP4
    etark1(i,j) = etark2(i,j) = a->eta(i,j);
    
    pgc->gcsl_start4(p,etark1,1);
    pgc->gcsl_start4(p,etark2,1);

}

void nhflow_fsf_fsm::step2(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, slice& etark1, slice &etark2, double alpha)
{

}

void nhflow_fsf_fsm::step3(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, slice& etark1, slice &etark2, double alpha)
{

}


