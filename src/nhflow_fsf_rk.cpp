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

#include"nhflow_fsf_rk.h"
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

nhflow_fsf_rk::nhflow_fsf_rk(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, patchBC_interface *ppBC) : epsi(p->A440*p->DXM), hx(p), hy(p)
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

void nhflow_fsf_rk::start(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow)
{
}

void nhflow_fsf_rk::ltimesave(lexer* p, fdm *a, slice &ls)
{
}

void nhflow_fsf_rk::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
    pupdate->start(p,a,pgc);
}

void nhflow_fsf_rk::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}

void nhflow_fsf_rk::step1(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, slice& etark1, slice &etark2, double alpha)
{

    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += u(i,j,k)*p->DZN[KP];

    VLOOP
	a->Q(i,j) += v(i,j,k)*p->DZN[KP];
    
    pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
    
    phxy->start(p,hx,hy,a->depth,a->wet,a->eta,a->P,a->Q);
    
    SLICELOOP1
    a->P(i,j) *= hx(i,j);

    SLICELOOP2
	a->Q(i,j) *= hy(i,j);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    SLICELOOP4
    etark1(i,j) =      a->eta(i,j)

                -      p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,etark1);
    pgc->gcsl_start4(p,etark1,1);
    
    p->sigma_update(p,a,pgc,etark1,a->eta,1.0);
    //p->omega_update(p,a,pgc,u,v,w,etark1,a->eta,1.0);
}

void nhflow_fsf_rk::step2(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, slice& etark1, slice &etark2, double alpha)
{
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += u(i,j,k)*p->DZN[KP];

    VLOOP
	a->Q(i,j) += v(i,j,k)*p->DZN[KP];
    
    pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
    
    phxy->start(p,hx,hy,a->depth,a->wet,etark1,a->P,a->Q);
    
    SLICELOOP1
    a->P(i,j) *= hx(i,j);

    SLICELOOP2
	a->Q(i,j) *= hy(i,j);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    SLICELOOP4
    etark2(i,j) = 0.75*a->eta(i,j) + 0.25*etark1(i,j)

                - 0.25*p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,etark2);
    pgc->gcsl_start4(p,etark2,1);
    
    p->sigma_update(p,a,pgc,etark2,etark1,0.25);
    //p->omega_update(p,a,pgc,u,v,w,etark2,etark1,0.25);
}

void nhflow_fsf_rk::step3(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, field &u, field&v, field&w, slice& etark1, slice &etark2, double alpha)
{
    // fill eta_n
    SLICELOOP4
    {
    a->eta_n(i,j) = a->eta(i,j);
    a->WL_n(i,j) = a->WL(i,j);
    }
    
    pgc->gcsl_start4(p,a->eta_n,1);
    
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

    ULOOP
    a->P(i,j) += u(i,j,k)*p->DZN[KP];

    VLOOP
	a->Q(i,j) += v(i,j,k)*p->DZN[KP];
    
    pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
    
    phxy->start(p,hx,hy,a->depth,a->wet,etark2,a->P,a->Q);
    
    SLICELOOP1
    a->P(i,j) *= hx(i,j);

    SLICELOOP2
	a->Q(i,j) *= hy(i,j);

	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);

    SLICELOOP4
    a->eta(i,j) = (1.0/3.0)*a->eta(i,j) + (2.0/3.0)*etark2(i,j)

                - (2.0/3.0)*p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);

    pflow->eta_relax(p,pgc,a->eta);
    pgc->gcsl_start4(p,a->eta,1);
    
    SLICELOOP4
    a->WL(i,j) = MAX(0.0, a->eta(i,j) + p->wd - a->bed(i,j));
    
    p->sigma_update(p,a,pgc,a->eta,etark2,2.0/3.0);
    //p->omega_update(p,a,pgc,u,v,w,a->eta,etark2,2.0/3.0);
}


