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

fnpf_breaking::fnpf_breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc) : bx(p), by(p), eta_t(p), t_break(p), B_coeff(p)
{
    // Default parameter values
    eta_I    = 0.65;    // onset threshold coefficient
    eta_F    = 0.15;    // cessation threshold coefficient
    delta_b  = p->A365;     // mixing length coefficient
    T_star   = 5.0;     // ramp-up time coefficient
    alpha_eta = 0.0;    // KFSBC dissipation (0 = off)


    ini_done = 0;
}

fnpf_breaking::~fnpf_breaking()
{
}

void fnpf_breaking::breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    if(p->A350==1)
    breaking_f(p, c, pgc, eta, eta_n, Fifsf, alpha);
    
    if(p->A350==2)
    breaking_kennedy(p, c, pgc, eta, eta_n, Fifsf, alpha);
}

