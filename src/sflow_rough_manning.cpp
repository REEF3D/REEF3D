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

#include"sflow_rough_manning.h"
#include"lexer.h"
#include"fdm2D.h"

#define WLIJ (WL(i,j)>p->A244?WL(i,j):1.0e20)

sflow_rough_manning::sflow_rough_manning(lexer* p) 
{
}

sflow_rough_manning::~sflow_rough_manning()
{
}

// Manning friction, tau_b/rho = g n^2 |u| u / h^(1/3), n = ks^(1/6)/20
// implicit in the momentum of the stage: d(UH)/dt = -gamma UH, gamma = g n^2 |u| / h^(4/3) from the stage state.
// The stage UHo = a UHn + (1-a)(UHs + dt L(UHs)) (SSP-RK) takes its Euler step as (UHs + dt L)/(1 + dt gamma);
// factor returns fr = 1/(1 + dt gamma).
// As an explicit source the friction is unstable for dt gamma > 2: thin water (gamma ~ h^(-4/3)) with the time step
// of the deep cells, e.g. 2 cm of sheet flow on a floodplain with 20 m cells (dt gamma 2.8: velocities of +-2 m/s
// instead of 0.05 m/s) and the film of a wetting front. The steady state is unchanged (L = gamma UH).

void sflow_rough_manning::factor(lexer *p, fdm2D *b, slice &fr, slice &U, slice &V, slice &WL)
{
    SLICELOOP4
    {
    fr(i,j) = 1.0;

        WETDRY
        {
        manning = pow(b->ks(i,j),1.0/6.0)/20.0;
        cf = pow(manning,2.0)*fabs(p->W22)/pow(WLIJ,1.0/3.0);

        fr(i,j) = 1.0/(1.0 + p->dt*cf*sqrt(U(i,j)*U(i,j) + V(i,j)*V(i,j))/WLIJ);
        }
    }
}
