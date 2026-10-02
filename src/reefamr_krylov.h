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
//
//  Restart: a space that also implements
//
//    void op_restart()                 RH = R, PV = 0, VV = 0 (leaf cells)
//
//  is restarted when the residual has not dropped by 5 % over RESTART_WIN iterations.
//  BiCGStab can stagnate with a constant residual (FNPF AMR: the roll and pitch unit-mode psi
//  solves of a floating body stalled at 1e-8..3e-7 relative for N 46 iterations); a new shadow
//  residual RH = R often recovers it within a few iterations (restarts receives their number).
//  If it stagnates again after a restart and the relative residual is below stagtol (0: never),
//  the solve stops there (stalls receives their number) instead of running to maxiter.

// default vector indices (a module may use its own)
enum { REEFAMR_NR=0, REEFAMR_NRH, REEFAMR_NPV, REEFAMR_NVV, REEFAMR_NS, REEFAMR_NT, REEFAMR_NPH, REEFAMR_NSH, REEFAMR_NTMP, REEFAMR_NVEC };

template<class S>
int reefamr_bicgstab(S &sp, double tol, int maxiter, double &bn, double &rn, double *tprec=nullptr, double *tapply=nullptr,
                     int *restarts=nullptr, double stagtol=0.0, int *stalls=nullptr)
{
    const int RESTART_WIN = 8;

    sp.op_start();

    bn = sqrt(sp.dot(S::B,S::B));
    rn = sqrt(sp.dot(S::R,S::R));
    double rho=1.0, alp=1.0, om=1.0;
    int it=0;
    double rbest=rn;
    int ibest=0, nostep=0;

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

        if constexpr (requires { sp.op_restart(); })
        {
            if(rn<0.95*rbest)
            {
                rbest = rn;
                ibest = it;
                nostep = 0;
            }
            else if(it-ibest>=RESTART_WIN && rn/bn>tol)
            {
                // stagnated again after a restart, close enough: stop
                if(nostep>0 && rn/bn<=stagtol)
                {
                    if(stalls)
                    ++(*stalls);
                    break;
                }

                ++nostep;
                sp.op_restart();
                rho = alp = om = 1.0;
                rbest = rn;
                ibest = it;
                if(restarts)
                ++(*restarts);
            }
        }
    }

    return it;
}

#endif
