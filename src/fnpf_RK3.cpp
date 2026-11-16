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

#include"fnpf_RK3.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fnpf_laplace_cds2.h"
#include"fnpf_fsf.h"
#include"fnpf_body.h"
#include"fnpf_amr.h"

fnpf_RK3::fnpf_RK3(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_rk(p,c,pgc),
                                                      erk1(p),erk2(p),frk1(p),frk2(p),ek(p),fk(p)
{
}

fnpf_RK3::~fnpf_RK3()
{
}

void fnpf_RK3::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini)
{
    // A 430: the tolerance of the stage solves (fnpf_laplace_cds2)
    fnpf_laplace_cds2::stage_tol = true;

    // A 321: start values of the stage solves extrapolated in time (1, 2), the coarse levels of REEFMG of the
    // first stage solve kept in the other two (2). Not with mesh refinement or resolved bodies, whose solves
    // are not those of plap0 alone.
    const bool sv = (p->A321>=1 && pamr==nullptr && p->X10==0);
    const bool kc = (p->A321>=2 && pamr==nullptr && p->X10==0 && p->N10==1);

    SLICELOOP4
    en(i,j) = c->eta(i,j);

// Step 1: t^n + dt
    stage_tendency(p,c,pgc,psolv,c->eta,ek,fk,0);

    SLICELOOP4
	erk1(i,j) = c->eta(i,j) + p->dt*ek(i,j);

    SLICELOOP4
	frk1(i,j) = c->Fifsf(i,j) + p->dt*fk(i,j);

    stage_surface(p,c,pgc,pflow,erk1,frk1,1.0,p->simtime+1.0*p->dt,1.0,1.0,0,true);
    stage_solve(p,c,pgc,psolv,erk1,frk1,sv?1:0,1.0,false);     // start value: Fi at t^n, for t^n + dt

// Step 2: t^n + dt/2
    stage_tendency(p,c,pgc,psolv,erk1,ek,fk,1);

    SLICELOOP4
	erk2(i,j) = 0.75*c->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*ek(i,j);

    SLICELOOP4
	frk2(i,j) = 0.75*c->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*fk(i,j);

    stage_surface(p,c,pgc,pflow,erk2,frk2,0.25,p->simtime+0.5*p->dt,1.0,0.5,1,false);
    stage_solve(p,c,pgc,psolv,erk2,frk2,sv?3:0,0.0,kc);        // start value: Fi at t^n + dt, for t^n + dt/2

// Step 3: t^n+1
    stage_tendency(p,c,pgc,psolv,erk2,ek,fk,2);

    SLICELOOP4
	c->eta(i,j) = (1.0/3.0)*c->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*ek(i,j);

    SLICELOOP4
	c->Fifsf(i,j) = (1.0/3.0)*c->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*fk(i,j);

    stage_surface(p,c,pgc,pflow,c->eta,c->Fifsf,2.0/3.0,p->simtime+1.0*p->dt,1.0,1.0,2,false);
    stage_solve(p,c,pgc,psolv,c->eta,c->Fifsf,sv?2:0,0.0,kc);  // start value: Fi at t^n + dt/2, for t^n + dt

    step_velocity(p,c,pgc);

    fnpf_laplace_cds2::stage_tol = false;

    if(pamr!=nullptr)
    pamr->step_end(p,c,pgc);
}

// mesh refinement (G 1): the Laplace equation is solved on all grids together
void fnpf_RK3::attach_amr(fnpf_amr *a)
{
    pamr = a;
    a->attach_level0(pf,this,this,this);
    a->attach_body(pbody);
    pbody->amr_attach(a);
    plap = a->laplace(plap,plap0);
}

void fnpf_RK3::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini
}
