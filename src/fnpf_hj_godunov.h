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

#ifndef FNPF_HJ_GODUNOV_H_
#define FNPF_HJ_GODUNOV_H_

// Godunov selection of the eta gradient in the kinematic FSBC (A315 3).
//
//   eta_t = Fz + sum_d h_d(eta_d),   h_d(q) = -F_d*q + Fz*q^2
//
// is a Hamilton-Jacobi equation in eta; written as eta_t + G(eta_x) = 0 with
// G(q) = F*q - Fz*q^2 per direction (separable).  With the left- and right-biased
// gradients qm, qp the Godunov Hamiltonian is
//
//   qm <= qp :  min of G over [qm,qp]        qm > qp :  max of G over [qp,qm]
//
// and G(q_sel) = G_godunov for the q_sel returned here, so the usual expression
// -F*q + Fz*(1 + q^2) evaluated with q_sel gives the Godunov flux.
// The upwinding by the characteristic speed G'(q) = F - 2*Fz*q of the previous
// stage (A315 1) picks one of qm, qp and misses the sonic point of the parabola:
// at a grid-scale dip with Fz<0 (or a spike with Fz>0) it takes the steeper side,
// Fz*q^2 then deepens the dip, the dip steepens, and the column runs away
// (moored box, N 61 stop).  Godunov gives q=0 at the bottom of the dip: the dip
// sinks with Fz only, as the exact (viscosity) solution.  Away from sonic points
// and for small slopes it is identical to A315 1.

inline double fnpf_hj_godunov(double qm, double qp, double F, double Fz)
{
    auto G = [&](double q) {return F*q - Fz*q*q;};
    
    double lo, hi;
    const bool minimize = (qm<=qp);
    
    if(minimize)
    {
    lo = qm;
    hi = qp;
    }
    else
    {
    lo = qp;
    hi = qm;
    }
    
    double qs = (G(lo)<=G(hi)) == minimize ? lo : hi;
    
    // vertex of the parabola: minimum for Fz<0, maximum for Fz>0
    if((minimize && Fz<0.0) || (!minimize && Fz>0.0))
    {
        const double qv = F/(2.0*Fz);
        
        if(qv>lo && qv<hi)
        qs = qv;
    }
    
    return qs;
}

#endif
