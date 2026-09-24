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

Discretisation: first-order upwind along s_hat,

        qb_P = ( qbe_P + |cx|*qb_upx + |cy|*qb_upy )/( 1 + |cx| + |cy| )

with cx = Ls*sgx/dx_upwind. The system matrix is a strictly diagonally
dominant M-matrix, so every update is a convex combination of qbe and
the upwind neighbours: unconditionally stable for any Ls/dx, maximum
principle, no negative qb.

Solver: the equation is steady, so it is solved to convergence in every
sediment step. Gauss-Seidel with alternating sweep directions (fast
sweeping): whenever the sweep direction matches the transport direction
the whole subdomain is solved in a single pass; the four orderings cover
all quadrants of s_hat. Ghostcell exchange after every sweep, iteration
stops once the global max update is below tol*max(qbe). The previous
Jacobi version stopped after a fixed number of sweeps, which turned the
steady relaxation into a pseudo-time lag that depended on the number of
sediment steps (S44, dtsed, RK stages) instead of on the physics.

Where the upwind neighbour is not a sediment cell (structure, inflow
ghost cell without bed) its term is dropped. This is the same fixed
point as a zero-gradient neighbour, but converges faster.

The relaxation is carried out on the bedload-only field qbn. s->qb is
only written at the end, because susp_qs() later adds the suspended load
qbs to s->qb, which must not re-enter the bedload relaxation.

Control:
  S33   0: equilibrium
        1: non-equilibrium, constant adaptation length Ls = S40 [m]
        2: non-equilibrium, van Rijn (1984) saltation length
           Ls = 3 d50 D*^0.6 T^0.9
        3: non-equilibrium, Phillips & Sutherland (1989)
           Ls = 4000 (theta - theta_cr) d50
        4: non-equilibrium, depth scaled Ls = S40*h   (S40 [-])
  S40   adaptation length [m] for S33 = 1, factor [-] for S33 = 4
  S49   maximum number of Gauss-Seidel sweeps per sediment step

Note on S33 = 2/3: both formulas give a single saltation hop, i.e.
O(1-100) d50. On typical grids Ls/dx << 1 and the result is practically
identical to equilibrium transport.
--------------------------------------------------------------------*/

void sediment_exner::non_equillibrium_solve(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    const double d50 = p->S20;
    const double visc = p->W2;
    const double grav = 9.81;
    const double Rstar = (p->S22 - p->W1)/p->W1;
    const double Dstar = d50*pow(fabs(Rstar)*grav/(visc*visc),1.0/3.0);
    const double ydir = p->y_dir;
    const int itermax = MAX(p->S49,1);
    const double tol = 1.0e-6;

    double uvel,vvel,umag;
    double sgx,sgy;
    double ustar2,ucrit2,Ti;
    double Lsc,cx,cy,num,den,qnew;
    double qbe_max=0.0;
    double Ls_max=0.0;
    double dqmax=0.0;
    int iter=0;


    // coefficients: transport direction and adaptation length ------------
    SEDSLICELOOP
    {
        // transport direction
        uvel = 0.5*(s->P(i,j) + s->P(i-1,j));
        vvel = 0.5*(s->Q(i,j) + s->Q(i,j-1))*ydir;

        umag = sqrt(uvel*uvel + vvel*vvel);

        sgx = umag>1.0e-10?uvel/umag:0.0;
        sgy = umag>1.0e-10?vvel/umag:0.0;


        // adaptation length
        Lsc = 0.0;

        if(p->S33==1)
        Lsc = p->S40;

        if(p->S33==2)
        {
        ucrit2 = s->shearvel_crit(i,j)*s->shearvel_crit(i,j);
        ustar2 = s->shearvel_eff(i,j)*s->shearvel_eff(i,j);

        Ti = ucrit2>1.0e-20?MAX((ustar2-ucrit2)/ucrit2,0.0):0.0;

        Lsc = 3.0*d50*pow(Dstar,0.6)*pow(Ti,0.9);
        }

        if(p->S33==3)
        Lsc = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j),0.0)*d50;

        if(p->S33==4)
        Lsc = p->S40*MAX(s->waterlevel(i,j),0.0);

        Lsc = MAX(Lsc,0.0);

        Ls_max = MAX(Ls_max,Lsc);


        // signed upwind coefficients, zero where the upwind cell carries no bed
        cx = 0.0;
        cy = 0.0;

        if(sgx>0.0 && p->DFBED[Im1J]>0 && p->DXP[IM1]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IM1];

        if(sgx<0.0 && p->DFBED[Ip1J]>0 && p->DXP[IP]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IP];

        if(sgy>0.0 && p->DFBED[IJm1]>0 && p->DYP[JM1]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JM1];

        if(sgy<0.0 && p->DFBED[IJp1]>0 && p->DYP[JP]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JP];

        cxn(i,j) = cx;
        cyn(i,j) = cy;

        qbe_max = MAX(qbe_max,s->qbe(i,j));
    }

    qbe_max = pgc->globalmax(qbe_max);
    Ls_max = pgc->globalmax(Ls_max);


    // initial guess: equilibrium, afterwards the previous sediment step
    if(noneq_ini==0)
    {
    SEDSLICELOOP
    qbn(i,j) = s->qbe(i,j);

    noneq_ini=1;
    }

    pgc->gcsl_start4(p,qbn,1);


    // Gauss-Seidel, alternating sweep directions ----------------------------
    for(iter=0; iter<itermax; ++iter)
    {
        const int sweep = iter%4;
        const bool xrev = (sweep==1 || sweep==2);
        const bool yrev = (sweep==2 || sweep==3);

        dqmax = 0.0;

        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        {
            i = xrev?(p->knox-1-ii):ii;
            j = yrev?(p->knoy-1-jj):jj;

            if(p->flagslice4[IJ]>0 && p->DFBED[IJ]>0)
            {
            cx = cxn(i,j);
            cy = cyn(i,j);

            num = s->qbe(i,j);
            den = 1.0;

            if(cx>0.0)
            {
            num += cx*qbn(i-1,j);
            den += cx;
            }

            if(cx<0.0)
            {
            num -= cx*qbn(i+1,j);
            den -= cx;
            }

            if(cy>0.0)
            {
            num += cy*qbn(i,j-1);
            den += cy;
            }

            if(cy<0.0)
            {
            num -= cy*qbn(i,j+1);
            den -= cy;
            }

            qnew = num/den;

            dqmax = MAX(dqmax,fabs(qnew-qbn(i,j)));

            qbn(i,j) = qnew;
            }
        }

        pgc->gcsl_start4(p,qbn,1);

        dqmax = pgc->globalmax(dqmax);

        if(dqmax<=tol*qbe_max)
        {
        ++iter;
        break;
        }
    }


    SEDSLICELOOP
    s->qb(i,j) = qbn(i,j);

    if(p->mpirank==0)
    cout<<"non-eq. bedload  Ls_max: "<<setprecision(4)<<Ls_max<<"  sweeps: "<<iter<<"  dq/qbe_max: "<<setprecision(3)<<(qbe_max>1.0e-20?dqmax/qbe_max:0.0)<<endl;
}
