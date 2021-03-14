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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"heat.h"
#include"concentration.h"
#include"momentum.h"
#include"sflow_eta_weno.h"
#include"sflow_hxy_weno.h"

nhflow_fsf_f::nhflow_fsf_f(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow) : 
                epsi(p->A440*p->DXM),depth(p),bed(p),L(p),hp(p),hx(p),hy(p)
{
	//peta = new sflow_eta_weno(p);
	//phxy = new sflow_hxy_weno(p);

	pupdate = new fluid_update_void();

    //pupdate = new fluid_update_fsf(p,a,pgc);
}

nhflow_fsf_f::~nhflow_fsf_f()
{
}

void nhflow_fsf_f::start(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow)
{
    
    // Momentum

    
    // fill eta_n
    SLICELOOP4
	{
    a->eta_n(i,j) = a->eta(i,j);
	L(i,j)=0.0;
	}
    pgc->gcsl_start4(p,a->eta_n,gcval_phi);
    
    
    // Calculate Eta
    SLICELOOP1
    a->P(i,j)=0.0;
    
    SLICELOOP2
    a->Q(i,j)=0.0;

	
    IULOOP
	JULOOP
    {
		KULOOP
        UCHECK
		{
			a->P(i,j) += a->u(i,j,k)*p->DZN[KP]*p->sigz[IJ];
		}
    }
    
    IVLOOP
	JVLOOP
	{	
			KVLOOP
            VCHECK
			{
            a->Q(i,j) += a->v(i,j,k)*p->DZN[KP]*p->sigz[IJ];
			}
    }
	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
	
	// depth update
	/*
    SLICELOOP4
	depth(i,j) = p->wd - bed(i,j);
	
	pgc->gcsl_start4(p,depth,50);
	
	SLICELOOP4
	hp(i,j) = a->eta(i,j) + p->wd - bed(i,j);
	
	pgc->gcsl_start4(p,hp,50);
	
	phxy->start(p,hx,hy,depth,a->eta,a->P,a->Q);
	
	pgc->gcsl_start1(p,hx,50);
	pgc->gcsl_start2(p,hy,50);
	
	SLICELOOP1
	a->P(i,j)/=hx(i,j);
	
	SLICELOOP2
	a->Q(i,j)/=hy(i,j);
	
	pgc->gcsl_start1(p,a->P,10);
    pgc->gcsl_start2(p,a->Q,11);
	
	// eta disc
	peta->start(p,a->eta,4,a->P,a->Q,depth,L);
    
	
	SLICELOOP4
	a->eta(i,j) +=	p->dt*L(i,j);	*/

    SLICELOOP4
    a->eta(i,j) -= p->dt*((a->P(i,j)-a->P(i-1,j))/p->DXN[IP] + (a->Q(i,j)-a->Q(i,j-1))/p->DYN[JP]);	  

    
    
    pflow->eta_relax(p,pgc,a->eta);
    
    pgc->gcsl_start4(p,a->eta,gcval_phi);
    
}

void nhflow_fsf_f::ltimesave(lexer* p, fdm *a, slice &ls)
{
}

void nhflow_fsf_f::update(lexer *p, fdm *a, ghostcell *pgc, slice &f)
{
    pupdate->start(p,a,pgc);
}

void nhflow_fsf_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}


