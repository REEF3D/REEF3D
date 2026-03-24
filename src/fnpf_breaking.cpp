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

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

fnpf_breaking::fnpf_breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc) : bx(p), by(p), eta_t(p), t_break(p), B_coeff(p),
                 S_local(p), b_strength(p), c_phase(p), k_local(p), L_break(p), B_ramp(p), eta_old(p), brk_flag(p),k_smooth(p)

{
    // Default parameter values
    eta_I    = 0.65;    // onset threshold coefficient
    eta_F    = 0.15;    // cessation threshold coefficient
    delta_b  = 1.2;     // mixing length coefficient
    T_star   = 5.0;     // ramp-up time coefficient
    alpha_eta = 0.0;    // KFSBC dissipation (0 = off)
    
    
    // Default parameters
    a_coeff       = 0.5;
    S0            = 0.08;
    alpha_L       = 0.35;
    T_star        = 3.0;
    alpha_eta     = 0.0;
    vb_max      = 10.0;
    eta_t_thresh  = 0.0;
    steep_method   = 2;
    cspeed_method  = 1;
    
    SLICELOOP4
    {
        c->vb(i,j)       = 0.0;
        S_local(i,j)    = 0.0;
        b_strength(i,j) = 0.0;
        c_phase(i,j)    = 1.0;
        k_local(i,j)    = 0.01;
        L_break(i,j)    = 1.0;
        B_ramp(i,j)     = 0.0;
        t_break(i,j)    = -1.0e20;
        eta_old(i,j)    = c->eta(i,j);
        eta_t(i,j)      = 0.0;
        brk_flag(i,j)   = 0;
    }

    pgc->gcsl_start4(p, c->vb, 1);
    pgc->gcsl_start4(p, S_local, 1);
    pgc->gcsl_start4(p, c_phase, 1);
    pgc->gcsl_start4(p, k_local, 1);


    ini_done = 0;
}

fnpf_breaking::~fnpf_breaking()
{
}

void fnpf_breaking::breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    if(p->A350==1 && p->A343==0)
    breaking_baquet_wd(p, c, pgc, eta, eta_n, Fifsf, alpha);
    
    if(p->A350==1 && p->A343==1)
    breaking_baquet_wd(p, c, pgc, eta, eta_n, Fifsf, alpha);
    
    if(p->A350==2)
    breaking_kennedy(p, c, pgc, eta, eta_n, Fifsf, alpha);
    
    if(p->A350==3)
    breaking_romero(p, c, pgc, eta, eta_n, Fifsf, alpha);
}

