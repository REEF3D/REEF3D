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

#ifndef FNPF_MEMORY_H_
#define FNPF_MEMORY_H_

#include"lexer.h"
#include"reefcouple.h"

// Which 3D fields of FNPF are stored. The core of the time step (Laplace assembly, velocities, bed
// BC, time step, free-surface output) needs none of these as 3D arrays; they are kept for the
// readers named below. (About 1.5 GB per array and 1e8 cells.)

// sig, sigx, sigy, sigxx: force gauges P 85, kinematics P 88, runup P 140, mesh refinement G 1,
// ice A 380 (fnpf_sigma_col.h)
inline bool fnpf_sigma_3D(lexer *p)
{
    return p->P85>0 || p->P88>0 || p->P140>0 || p->G1>0 || p->A380>0;
}

// U, V, W: floating bodies X 10, force gauges P 85, kinematics P 88, runup P 140, state files P 40,
// mesh refinement G 1, ice A 380, CORAL Z 40, checkpoints P 48 / I 45, the runtime coupling to a
// receiver; otherwise the surface values Us, Vs, Ws, and the 3D field is made for a volume output
// (fdm_fnpf::vel3D_make)
inline bool fnpf_vel_3D(lexer *p)
{
    return p->X10>0 || p->P85>0 || p->P88>0 || p->P140>0 || p->P40>0 || p->G1>0 || p->A380>0
        || p->Z40>0 || p->P48>0 || p->P48_dt>0.0 || p->I45==1 || reefcouple::source();
}

// Uin: Dirichlet wave generation B 98 3, 4, the active beaches B 99 > 2, floating bodies X 10,
// mesh refinement G 1, checkpoints; elsewhere it stays zero (gc_fivec reads it as zero)
inline bool fnpf_uin_3D(lexer *p)
{
    return p->B98>2 || p->B99>2 || p->X10>0 || p->G1>0 || p->P48>0 || p->P48_dt>0.0 || p->I45==1;
}

#endif
