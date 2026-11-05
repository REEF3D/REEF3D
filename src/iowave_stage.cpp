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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"iowave.h"
#include"lexer.h"
#include"fdm2D.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>

//  B 95: step-size consistent relaxation zones (SFLOW and NHFLOW, B 98 2).
//
//  The relaxation zones act in every RK stage: U = (1-r) U_target + r U.  Without B 95 the
//  targets are taken once per step at its start and r does not depend on dt, so the generated
//  wave depends on the number of steps per wave period (and lags the stage outputs by up to dt).
//  With B 95 the momentum stage calls wavegen_*_stage before the stage: the targets at the time
//  of the stage output (RK3: t+dt, t+dt/2, t+dt; RK2: t+dt, t+dt) and the relaxation factors
//  r^(dt/dt_ref) of the generation and beach zones, dt_ref = B 95 (> 0) or, with B 95 -1, the
//  step of the finest level (dt/2^G1 with mesh refinement subcycling G 7 1, else dt).  Over a
//  fixed time the zones then relax by the same amount whatever the step.

void iowave::relax_stage_factor(lexer *p)
{
    double dtref = p->B95;
    if(p->B95<0.0)
    dtref = (p->G7==1 && p->G1>0) ? p->dt/double(1<<p->G1) : p->dt;

    const double f = p->dt/(dtref>1.0e-20 ? dtref : 1.0e-20);

    if(relax4_wg0==nullptr)
    {
        relax4_wg0 = new slice4(p);
        relax4_nb0 = new slice4(p);
        SLICEBASELOOP
        {
        (*relax4_wg0)(i,j) = relax4_wg(i,j);
        (*relax4_nb0)(i,j) = relax4_nb(i,j);
        }
        relax_fac = 1.0;
    }

    if(f==relax_fac)
    return;

    relax_fac = f;
    SLICEBASELOOP
    {
    const double rw = (*relax4_wg0)(i,j), rb = (*relax4_nb0)(i,j);
    relax4_wg(i,j) = (rw>0.0 && rw<1.0) ? pow(rw,f) : rw;
    relax4_nb(i,j) = (rb>0.0 && rb<1.0) ? pow(rb,f) : rb;
    }
}

void iowave::wavegen_stage_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double t)
{
    if(p->B95==0.0 || p->B98!=2)
    return;

    relax_stage_factor(p);

    // the relaxation precalc of wavegen_precalc_nhflow at time t (it takes the time from simtime)
    const double ts = p->simtime;
    p->simtime = t;

    starttime=pgc->timer();

    if(p->B89==0)
    nhflow_precalc_relax(p,d,pgc);

    if(p->B89==1)
    {
    nhflow_wavegen_precalc_decomp_time(p,pgc);
    nhflow_wavegen_precalc_decomp_relax(p,d,pgc);
    }

    p->wavecalctime+=pgc->timer()-starttime;

    p->wavetime = t;
    p->simtime = ts;
}

void iowave::wavegen_2D_stage(lexer *p, fdm2D *b, ghostcell *pgc, double t)
{
    if(p->B95==0.0 || p->B98!=2)
    return;

    relax_stage_factor(p);

    const double ts = p->simtime;
    p->simtime = t;
    wavegen_2D_precalc(p,b,pgc);
    p->wavetime = t;
    p->simtime = ts;
}
