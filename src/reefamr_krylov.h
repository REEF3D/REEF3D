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

#ifndef REEFAMR_KRYLOV_H_
#define REEFAMR_KRYLOV_H_

#include<cmath>
#include<mpi.h>

//  BiCGStab, right preconditioned, for a composite elliptic problem on the leaf cells of
//  a REEFAMR hierarchy.  The module supplies the vector space S: the vectors live on all
//  grids of the hierarchy and are named by the indices below (-1: the solution, which holds
//  the initial guess), and S implements
//
//    void apply(int x, int y)          y = A x on the leaf cells
//    void prec(int r, int z)           z = M^-1 r (e.g. one FAC sweep)
//    double dot(int a, int b)          global dot product over the leaf cells
//    void op_start()                   R = S - V (V = A x0, S = rhs), RH = R, P = 0, V = 0
//    void op_p(double beta, double om) P = R + beta*(P - om*V)
//    void op_s(double alp)             S = R - alp*V
//    void op_x(double alp, double om)  x += alp*PH + om*SH, R = S - om*T (leaf cells)
//
//  The right-hand side is expected in S (NS) and A x0 in V (NVV) before op_start.
//  Returns the iteration count; bn, rn receive the norms of the right-hand side and of the
//  final residual.  tprec, tapply accumulate the time in the preconditioner and the operator.

enum { REEFAMR_NR=0, REEFAMR_NRH, REEFAMR_NPV, REEFAMR_NVV, REEFAMR_NS, REEFAMR_NT, REEFAMR_NPH, REEFAMR_NSH, REEFAMR_NTMP, REEFAMR_NVEC };

template<class S>
int reefamr_bicgstab(S &sp, double tol, int maxiter, double &bn, double &rn, double *tprec=nullptr, double *tapply=nullptr)
{
    sp.op_start();

    bn = sqrt(sp.dot(REEFAMR_NS,REEFAMR_NS));
    rn = sqrt(sp.dot(REEFAMR_NR,REEFAMR_NR));
    double rho=1.0, alp=1.0, om=1.0;
    int it=0;

    if(bn>0.0)
    while(rn/bn>tol && it<maxiter)
    {
        double rhon = sp.dot(REEFAMR_NRH,REEFAMR_NR);
        double beta = (rhon/rho)*(alp/om);
        rho = rhon;

        sp.op_p(beta,om);

        double ta = MPI_Wtime();
        sp.prec(REEFAMR_NPV,REEFAMR_NPH);
        double tb = MPI_Wtime();
        sp.apply(REEFAMR_NPH,REEFAMR_NVV);
        if(tprec)  *tprec += tb-ta;
        if(tapply) *tapply += MPI_Wtime()-tb;
        alp = rho/sp.dot(REEFAMR_NRH,REEFAMR_NVV);

        sp.op_s(alp);

        ta = MPI_Wtime();
        sp.prec(REEFAMR_NS,REEFAMR_NSH);
        tb = MPI_Wtime();
        sp.apply(REEFAMR_NSH,REEFAMR_NT);
        if(tprec)  *tprec += tb-ta;
        if(tapply) *tapply += MPI_Wtime()-tb;
        double tt = sp.dot(REEFAMR_NT,REEFAMR_NT);
        om = tt>0.0 ? sp.dot(REEFAMR_NT,REEFAMR_NS)/tt : 0.0;

        sp.op_x(alp,om);

        rn = sqrt(sp.dot(REEFAMR_NR,REEFAMR_NR));
        ++it;
    }

    return it;
}

#endif
