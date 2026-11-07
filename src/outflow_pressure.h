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

#ifndef OUTFLOW_PRESSURE_H_
#define OUTFLOW_PRESSURE_H_

#include<cmath>
#include"lexer.h"
#include"fdm.h"
#include"heaviside.h"
#include"increment.h"

// pressure in the outflow ghost cells of the boundary cell (i,j,k), B 77
// (one definition for the ghost cells of ioflow_f / iowave / ioflow_v and the Poisson rows)
//
//  0  zero gradient: the cell value
//  1  controlled outflow, pressure condition: hydrostatic from the outflow water level
//     (F 60 / F 62) for F 50 = 2, 3; the cell value for F 50 = 1, 4 and with a relaxation
//     beach (B 99 = 1, 2), where the beach sets the free surface and the velocities
//  2  controlled outflow, free surface condition (F 62, ioflow_waterlevel): hydrostatic from the
//     outflow water level as well (F 62, ramped from F 60 over F 64), so that the pressure agrees
//     with the level set held at that level in the ghost cells; the cell value without F 62 and
//     with a relaxation beach (B 99 = 1, 2)
//
// p->fsfout: F 60 / F 62 (ini_phi), during the run the ramped F 62 (ioflow_f::fsfinflow).
// In the Poisson row of an outflow face (poisson_pcorr) the Dirichlet value of the correction is
// outflow_pressure() - cell pressure: the boundary face velocity is not part of the momentum step,
// so the ghost pressure alone has no effect.
// 10  free stream: zero pressure in the liquid, the cell value in the gas, (1-H) p
inline double outflow_pressure(lexer *p, fdm *a, int i, int j, int k)
{
    const int marge = increment::marge;

    const bool beach = (p->B90>0 && (p->B99==1 || p->B99==2));

    if(((p->B77==1 && (p->F50==2 || p->F50==3)) || (p->B77==2 && p->F62>-1.0e20)) && !beach)
    {
    double z = (p->G502==1) ? p->ZSP[IJK] : p->ZP[KP];

    return (p->fsfout - z)*a->ro(i,j,k)*fabs(p->W22);
    }

    if(p->B77==10)
    {
    double eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);

    return (1.0-heaviside(a->phi(i,j,k),eps))*a->press(i,j,k);
    }

    return a->press(i,j,k);
}

#endif
