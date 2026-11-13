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

#include"fnpf_RK4.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_laplace_cds2.h"
#include"fnpf_fsf.h"

fnpf_RK4::fnpf_RK4(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_rk(p,c,pgc),
                                                      erk1(p),erk2(p),erk3(p),erk4(p),
                                                      frk1(p),frk2(p),frk3(p),frk4(p),
                                                      erk(p),frk(p)
{
}

fnpf_RK4::~fnpf_RK4()
{
}

void fnpf_RK4::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini)
{
    // A 430: the tolerance of the stage solves (fnpf_laplace_cds2)
    fnpf_laplace_cds2::stage_tol = true;

    // A 321: start values of the stage solves extrapolated in time (1, 2), the coarse levels of REEFMG of the
    // first stage solve kept in the other three (2). Not with resolved bodies, whose solves are not
    // those of plap0 alone. The second and the fourth solve start from the solution of the solve
    // before them, which is at their time.
    const bool sv = (p->A321>=1 && p->X10==0);
    const bool kc = (p->A321>=2 && p->X10==0 && p->N10==1);

    // the relaxation zones of the stages with the weight of the stage: as strong per time as in RK3 (iowave_stage.cpp)

    SLICELOOP4
    en(i,j) = c->eta(i,j);

// Step 1: t^n + dt/2
    stage_tendency(p,c,pgc,psolv,c->eta,erk1,frk1,0);

    SLICELOOP4
    {
    erk(i,j) = c->eta(i,j)   + 0.5*p->dt*erk1(i,j);
    frk(i,j) = c->Fifsf(i,j) + 0.5*p->dt*frk1(i,j);
    }

    stage_surface(p,c,pgc,pflow,erk,frk,0.5,p->simtime+0.5*p->dt,11.0/12.0,0.5,0,true);
    stage_solve(p,c,pgc,psolv,erk,frk,sv?1:0,0.5,false);      // start value: Fi at t^n, for t^n + dt/2

// Step 2: t^n + dt/2
    stage_tendency(p,c,pgc,psolv,erk,erk2,frk2,1);

    SLICELOOP4
    {
    erk(i,j) = c->eta(i,j)   + 0.5*p->dt*erk2(i,j);
    frk(i,j) = c->Fifsf(i,j) + 0.5*p->dt*frk2(i,j);
    }

    stage_surface(p,c,pgc,pflow,erk,frk,0.5,p->simtime+0.5*p->dt,11.0/12.0,0.5,1,false);
    stage_solve(p,c,pgc,psolv,erk,frk,0,0.0,kc);

// Step 3: t^n + dt
    stage_tendency(p,c,pgc,psolv,erk,erk3,frk3,2);

    SLICELOOP4
    {
    erk(i,j) = c->eta(i,j)   + p->dt*erk3(i,j);
    frk(i,j) = c->Fifsf(i,j) + p->dt*frk3(i,j);
    }

    stage_surface(p,c,pgc,pflow,erk,frk,1.0,p->simtime+1.0*p->dt,11.0/6.0,1.0,2,false);
    stage_solve(p,c,pgc,psolv,erk,frk,sv?2:0,0.0,kc);         // start value: Fi at t^n + dt/2, for t^n + dt

// Step 4: t^n+1
    stage_tendency(p,c,pgc,psolv,erk,erk4,frk4,3);

    SLICELOOP4
    {
    c->eta(i,j)   = c->eta(i,j)   + p->dt*(1.0/6.0)*(erk1(i,j) + 2.0*erk2(i,j) + 2.0*erk3(i,j) + erk4(i,j));
    c->Fifsf(i,j) = c->Fifsf(i,j) + p->dt*(1.0/6.0)*(frk1(i,j) + 2.0*frk2(i,j) + 2.0*frk3(i,j) + frk4(i,j));
    }

    // the relaxation that enters the state, with the weight of the three stages of RK3
    stage_surface(p,c,pgc,pflow,c->eta,c->Fifsf,1.0,p->simtime+p->dt,11.0/6.0,1.0,3,false);
    stage_solve(p,c,pgc,psolv,c->eta,c->Fifsf,0,0.0,kc);

    step_velocity(p,c,pgc);

    fnpf_laplace_cds2::stage_tol = false;
}

void fnpf_RK4::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   // coastline ini

    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
}
