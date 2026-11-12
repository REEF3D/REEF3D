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
#include"sediment_curv.h"
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

Where the upwind neighbour has no erodible bed (a structure in the bed,
DFBED < 0, or a solid column inside the domain, no slice cell) it
supplies no bed load: its term stays in the equation with qb = 0, so
that the transport behind a structure, an apron or a sill grows from
zero over Ls (dropping the term gave qb = qbe in the first cell, the
equilibrium erosion at the edge). A prescribed bed load at the inflow
(S 129) is the upwind value of the first cells (the inflow ghost cells
of qbn, qb_inflow): the transport adapts from the supply to the capacity
over Ls. Without it the inflow ghost cells are zero gradient, qb = qbe.

The relaxation is carried out on the bedload-only field qbn. s->qb is
only written at the end, because susp_qs() later adds the suspended load
qbs to s->qb, which must not re-enter the bedload relaxation.

Control:
  S33   0: equilibrium
        1: non-equilibrium, constant adaptation length Ls = S40 [m]
        2: non-equilibrium, S40 saltation lengths of van Rijn (1984)
           Ls = S40 * 3 d50 D*^0.6 T^0.9   (S40 >= 1, number of hops)
        3: non-equilibrium, Phillips & Sutherland (1989)
           Ls = 4000 (theta - theta_cr) d50
        4: non-equilibrium, depth scaled Ls = S40*h   (S40 [-])
  S40   adaptation length [m] for S33 = 1, number of saltation lengths [-]
        for S33 = 2, factor [-] for S33 = 4
  S49   maximum number of Gauss-Seidel sweeps per sediment step

Note on S33 = 2: one saltation length is 10-200 d50 (7-100 mm for
d50 = 0.5 mm, T = 1-20), below the cell size of most grids, where the
result equals equilibrium transport; S40 is the number of hops a grain
makes before it comes to rest. S33 = 3 gives 4000 (theta - theta_cr) d50,
0.04-1 m for the same sand.
--------------------------------------------------------------------*/

void sediment_exner::non_equillibrium_solve(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    // curvilinear grid (G 50): sediment_curv.h
    if(p->geo_curv==1)
    {
    non_equillibrium_solve_curv(p,pgc,s);
    return;
    }

    const double d50 = s->dk;
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

        Lsc = p->S40*3.0*d50*pow(Dstar,0.6)*pow(Ti,0.9);
        }

        if(p->S33==3)
        Lsc = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j),0.0)*d50;

        if(p->S33==4)
        Lsc = p->S40*MAX(s->waterlevel(i,j),0.0);

        Lsc = MAX(Lsc,0.0);

        Ls_max = MAX(Ls_max,Lsc);


        // signed upwind coefficients (an upwind cell without a bed supplies no bed load: qb = 0 in the sweep)
        cx = 0.0;
        cy = 0.0;

        if(sgx>0.0 && p->DXP[IM1]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IM1];

        if(sgx<0.0 && p->DXP[IP]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IP];

        if(sgy>0.0 && p->DYP[JM1]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JM1];

        if(sgy<0.0 && p->DYP[JP]>1.0e-20)
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

    // mesh refinement and multi-block grids: the cells around the grid and the covered cells hold their values
    // (lagfix, sediment_exner_lag.cpp); elsewhere this is the ghost exchange
    lag_hold(p,pgc,qbn,0);
    qb_inflow(p,s,frac_k,qbn);


    // a neighbour without an erodible bed supplies no bed load: a structure in the bed (DFBED < 0) or a solid
    // column inside the domain (no slice cell: NHFLOW and SFLOW solids, piers). The cells beyond the edges of the
    // domain keep their ghost values (zero gradient, or the supply of the inflow)
    auto supply = [&](int di, int dj)
    {
        const int nb = (i+di-p->imin)*p->jmax + (j+dj-p->jmin);
        const int gi = i+di+p->origin_i;
        const int gj = j+dj+p->origin_j;
        
        if(p->DFBED[nb]<0)
        return false;
        
        if(gi>=0 && gi<p->gknox && gj>=0 && gj<p->gknoy && p->flagslice4[nb]<0)
        return false;
        
        return true;
    };

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

            if(p->flagslice4[IJ]>0 && p->DFBED[IJ]>0 && lagfix(i,j)<=0.0)
            {
            cx = cxn(i,j);
            cy = cyn(i,j);

            num = s->qbe(i,j);
            den = 1.0;

            if(cx>0.0)
            {
            num += cx*(supply(-1,0) ? qbn(i-1,j) : 0.0);
            den += cx;
            }

            if(cx<0.0)
            {
            num -= cx*(supply(1,0) ? qbn(i+1,j) : 0.0);
            den -= cx;
            }

            if(cy>0.0)
            {
            num += cy*(supply(0,-1) ? qbn(i,j-1) : 0.0);
            den += cy;
            }

            if(cy<0.0)
            {
            num -= cy*(supply(0,1) ? qbn(i,j+1) : 0.0);
            den -= cy;
            }

            qnew = num/den;

            dqmax = MAX(dqmax,fabs(qnew-qbn(i,j)));

            qbn(i,j) = qnew;
            }
        }

        lag_hold(p,pgc,qbn,1);
        qb_inflow(p,s,frac_k,qbn);

        dqmax = pgc->globalmax(dqmax);

        if(iter==0)
        {
        lag_d[2] = dqmax;
        lag_f[2] = qbe_max;
        }

        if(dqmax<=tol*qbe_max)
        {
        ++iter;
        break;
        }
    }


    SEDSLICELOOP
    s->qb(i,j) = qbn(i,j);

    if(p->mpirank==0 && lag_print)
    cout<<"non-eq. bedload  Ls_max: "<<setprecision(4)<<Ls_max<<"  sweeps: "<<iter<<"  dq/qbe_max: "<<setprecision(3)<<(qbe_max>1.0e-20?dqmax/qbe_max:0.0)<<endl;
}

// ---- curvilinear grid (G 50): the functions above with the geometry of sediment_curv.h

void sediment_exner::non_equillibrium_solve_curv(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    const double d50 = s->dk;
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
        if(p->geo_curv==1)
        sediment_curv::cellvel(p,s,i,j,uvel,vvel);
        else
        {
        uvel = 0.5*(s->P(i,j) + s->P(i-1,j));
        vvel = 0.5*(s->Q(i,j) + s->Q(i,j-1))*ydir;
        }

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

        Lsc = p->S40*3.0*d50*pow(Dstar,0.6)*pow(Ti,0.9);
        }

        if(p->S33==3)
        Lsc = 4000.0*MAX(s->shields_eff(i,j)-s->shields_crit(i,j),0.0)*d50;

        if(p->S33==4)
        Lsc = p->S40*MAX(s->waterlevel(i,j),0.0);

        Lsc = MAX(Lsc,0.0);

        Ls_max = MAX(Ls_max,Lsc);


        // signed upwind coefficients (an upwind cell without a bed supplies no bed load: qb = 0 in the sweep)
        cx = 0.0;
        cy = 0.0;

        if(sgx>0.0 && p->DXP[IM1]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IM1];

        if(sgx<0.0 && p->DXP[IP]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IP];

        if(sgy>0.0 && p->DYP[JM1]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JM1];

        if(sgy<0.0 && p->DYP[JP]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JP];

        // curvilinear grid (G 50): Ls s.grad(qb) upwinded per index direction, the contravariant
        // components of s (per index, from the centre differences over two cells)
        if(p->geo_curv==1)
        {
        const double sxi = sediment_curv::contra(p,1,i,j,sgx,sgy);
        const double seta = sediment_curv::contra(p,2,i,j,sgx,sgy);
        cx = cy = 0.0;
        if(sxi>0.0)
        cx = Lsc*sxi;
        if(sxi<0.0)
        cx = Lsc*sxi;
        if(seta>0.0)
        cy = Lsc*seta;
        if(seta<0.0)
        cy = Lsc*seta;
        }
        
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

    // mesh refinement and multi-block grids: the cells around the grid and the covered cells hold their values
    // (lagfix, sediment_exner_lag.cpp); elsewhere this is the ghost exchange
    lag_hold(p,pgc,qbn,0);
    qb_inflow(p,s,frac_k,qbn);


    // a neighbour without an erodible bed supplies no bed load: a structure in the bed (DFBED < 0) or a solid
    // column inside the domain (no slice cell: NHFLOW and SFLOW solids, piers). The cells beyond the edges of the
    // domain keep their ghost values (zero gradient, or the supply of the inflow)
    auto supply = [&](int di, int dj)
    {
        const int nb = (i+di-p->imin)*p->jmax + (j+dj-p->jmin);
        const int gi = i+di+p->origin_i;
        const int gj = j+dj+p->origin_j;
        
        if(p->DFBED[nb]<0)
        return false;
        
        if(gi>=0 && gi<p->gknox && gj>=0 && gj<p->gknoy && p->flagslice4[nb]<0)
        return false;
        
        return true;
    };

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

            if(p->flagslice4[IJ]>0 && p->DFBED[IJ]>0 && lagfix(i,j)<=0.0)
            {
            cx = cxn(i,j);
            cy = cyn(i,j);

            num = s->qbe(i,j);
            den = 1.0;

            if(cx>0.0)
            {
            num += cx*(supply(-1,0) ? qbn(i-1,j) : 0.0);
            den += cx;
            }

            if(cx<0.0)
            {
            num -= cx*(supply(1,0) ? qbn(i+1,j) : 0.0);
            den -= cx;
            }

            if(cy>0.0)
            {
            num += cy*(supply(0,-1) ? qbn(i,j-1) : 0.0);
            den += cy;
            }

            if(cy<0.0)
            {
            num -= cy*(supply(0,1) ? qbn(i,j+1) : 0.0);
            den -= cy;
            }

            qnew = num/den;

            dqmax = MAX(dqmax,fabs(qnew-qbn(i,j)));

            qbn(i,j) = qnew;
            }
        }

        lag_hold(p,pgc,qbn,1);
        qb_inflow(p,s,frac_k,qbn);

        dqmax = pgc->globalmax(dqmax);

        if(iter==0)
        {
        lag_d[2] = dqmax;
        lag_f[2] = qbe_max;
        }

        if(dqmax<=tol*qbe_max)
        {
        ++iter;
        break;
        }
    }


    SEDSLICELOOP
    s->qb(i,j) = qbn(i,j);

    if(p->mpirank==0 && lag_print)
    cout<<"non-eq. bedload  Ls_max: "<<setprecision(4)<<Ls_max<<"  sweeps: "<<iter<<"  dq/qbe_max: "<<setprecision(3)<<(qbe_max>1.0e-20?dqmax/qbe_max:0.0)<<endl;
}
