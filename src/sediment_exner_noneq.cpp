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

#include"sediment_exner.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"solver2D.h"
#include"sediment_exnerdisc.h"
#include<math.h>

/*--------------------------------------------------------------------
Non-equilibrium (non-capacity) bedload transport.

The bedload flux does not adjust instantaneously to the local transport
capacity qbe: grains entrained at one location travel a finite distance
before they are deposited again. The transport rate therefore relaxes
towards the capacity over the adaptation length Ls along the transport
path

        qb + Ls*( s_hat . grad(qb) ) = qbe ,     s_hat = (u,v)/|U|

with Ls -> 0 recovering the equilibrium closure qb = qbe.

Discretisation: first-order upwind along s_hat, solved for the diagonal

        qb_P = ( qbe_P + cx*qb_upwind + cy*qb_upwind )/( 1 + cx + cy )

with cx = Ls*|sgx|/dx_upwind >= 0. Each update is a convex combination
of qbe and the upwind neighbours, so the sweep is unconditionally stable
for any Ls/dx, satisfies a maximum principle and cannot generate
negative qb. Jacobi sweeps propagate the information one cell per sweep
with a ghostcell exchange in between; the solution of the previous
sediment time step (qbn) is used as the initial guess.

The relaxation is carried out on the bedload-only field qbn. s->qb is
only written at the end, because susp_qs() later adds the suspended load
qbs to s->qb, which must not re-enter the bedload relaxation.

Control:
  S33   0: equilibrium
        1: non-equilibrium, constant adaptation length Ls = S40
        2: non-equilibrium, van Rijn saltation length
        3: non-equilibrium, Phillips & Sutherland
  S40   adaptation length [m], used for S33 = 1
  S49   number of relaxation sweeps per sediment time step
--------------------------------------------------------------------*/

void sediment_exner::non_equillibrium_solve(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    const double d50 = p->S20;
    const double visc = p->W2;
    const double grav = 9.81;
    const double Rstar = (p->S22 - p->W1)/p->W1;
    const double Dstar = d50*pow(fabs(Rstar)*grav/(visc*visc),1.0/3.0);
    const double ydir = p->y_dir;
    const int itermax = 10;

    double uvel,vvel,umag;
    double sgx,sgy;
    double dxu,dyu,qxu,qyu,cx,cy;
    double ustar2,ucrit2,Ti;
    double Ls;


    // initial guess: equilibrium
    if(noneq_ini==0)
    {
    SEDSLICELOOP
    qbn(i,j) = s->qbe(i,j);

    pgc->gcsl_start4(p,qbn,1);

    noneq_ini=1;
    }

    SEDSLICELOOP
    q0(i,j) = qbn(i,j);

    pgc->gcsl_start4(p,q0,1);


    for(int qn=0; qn<itermax; ++qn)
    {
    SEDSLICELOOP
    {
        // transport direction
        uvel = 0.5*(s->P(i,j) + s->P(i-1,j));
        vvel = 0.5*(s->Q(i,j) + s->Q(i,j-1))*ydir;

        umag = sqrt(uvel*uvel + vvel*vvel);

        sgx = umag>1.0e-10?uvel/umag:0.0;
        sgy = umag>1.0e-10?vvel/umag:0.0;


        // adaptation length
        ucrit2 = s->shearvel_crit(i,j)*s->shearvel_crit(i,j);
        ustar2 = s->shearvel_eff(i,j)*s->shearvel_eff(i,j);

        Ti = ucrit2>1.0e-20?MAX((ustar2-ucrit2)/ucrit2,0.0):0.0;

        Ls = 0.1;

        if(p->S33==2)
        Ls = 3.0*d50*pow(Dstar,0.6)*pow(Ti,0.9);

        if(p->S33==3)
        Ls = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j),0.0)*d50;

        Ls = MAX(Ls,0.0);


        // upwind neighbors along s_hat; fall back to the local value
        // where there is no sediment bed (zero gradient)
        if(sgx>=0.0)
        {
        dxu = p->DXP[IM1];
        qxu = p->DFBED[Im1J]>0?q0(i-1,j):q0(i,j);
        }

        if(sgx<0.0)
        {
        dxu = p->DXP[IP];
        qxu = p->DFBED[Ip1J]>0?q0(i+1,j):q0(i,j);
        }

        if(sgy>=0.0)
        {
        dyu = p->DYP[JM1];
        qyu = p->DFBED[IJm1]>0?q0(i,j-1):q0(i,j);
        }

        if(sgy<0.0)
        {
        dyu = p->DYP[JP];
        qyu = p->DFBED[IJp1]>0?q0(i,j+1):q0(i,j);
        }

        cx = dxu>1.0e-20?Ls*fabs(sgx)/dxu:0.0;
        cy = dyu>1.0e-20?Ls*fabs(sgy)/dyu:0.0;


        qbn(i,j) = (s->qbe(i,j) + cx*qxu + cy*qyu)/(1.0 + cx + cy);
    }

    SEDSLICELOOP
    q0(i,j) = qbn(i,j);

    pgc->gcsl_start4(p,q0,1);
    }


    SEDSLICELOOP
    s->qb(i,j) = qbn(i,j);
}
