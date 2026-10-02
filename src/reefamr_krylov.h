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
//  grids of the hierarchy and are named by indices that S defines as static constants
//  (B: right-hand side, R, RH, PV, VV, S, T, PH, SH: work vectors), and S implements
//
//    void apply(int x, int y)          y = A x on the leaf cells
//    void prec(int r, int z)           z = M^-1 r (e.g. one FAC sweep, line solves)
//    double dot(int a, int b)          global dot product over the leaf cells
//    void op_start()                   R = B - VV (VV = A x0), RH = R, PV = 0, VV = 0
//    void op_p(double beta, double om) PV = R + beta*(PV - om*VV)
//    void op_s(double alp)             S = R - alp*VV
//    void op_x(double alp, double om)  x += alp*PH + om*SH, R = S - om*T (leaf cells)
//
//  The right-hand side is expected in B and A x0 in VV before the call (B may be the same
//  vector as S: it is read before S is first written).  Returns the iteration count; bn, rn
//  receive the norms of the right-hand side and of the final residual.  tprec, tapply
//  accumulate the time in the preconditioner and the operator.

// default vector indices (a module may use its own)
enum { REEFAMR_NR=0, REEFAMR_NRH, REEFAMR_NPV, REEFAMR_NVV, REEFAMR_NS, REEFAMR_NT, REEFAMR_NPH, REEFAMR_NSH, REEFAMR_NTMP, REEFAMR_NVEC };

template<class S>
int reefamr_bicgstab(S &sp, double tol, int maxiter, double &bn, double &rn, double *tprec=nullptr, double *tapply=nullptr)
{
    sp.op_start();

    bn = sqrt(sp.dot(S::B,S::B));
    rn = sqrt(sp.dot(S::R,S::R));
    double rho=1.0, alp=1.0, om=1.0;
    int it=0;

    if(bn>0.0)
    while(rn/bn>tol && it<maxiter)
    {
        double rhon = sp.dot(S::RH,S::R);
        double beta = (rhon/rho)*(alp/om);
        rho = rhon;

        sp.op_p(beta,om);

        double ta = MPI_Wtime();
        sp.prec(S::PV,S::PH);
        double tb = MPI_Wtime();
        sp.apply(S::PH,S::VV);
        if(tprec)  *tprec += tb-ta;
        if(tapply) *tapply += MPI_Wtime()-tb;
        alp = rho/sp.dot(S::RH,S::VV);

        sp.op_s(alp);

        ta = MPI_Wtime();
        sp.prec(S::S,S::SH);
        tb = MPI_Wtime();
        sp.apply(S::SH,S::T);
        if(tprec)  *tprec += tb-ta;
        if(tapply) *tapply += MPI_Wtime()-tb;
        double tt = sp.dot(S::T,S::T);
        om = tt>0.0 ? sp.dot(S::T,S::S)/tt : 0.0;

        sp.op_x(alp,om);

        rn = sqrt(sp.dot(S::R,S::R));
        ++it;
    }

    return it;
}

#endif
