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

#ifndef NHFLOW_STRATIFICATION_H_
#define NHFLOW_STRATIFICATION_H_

class lexer;
class fdm_nhf;

// ---------------------------------------------------------------------------
// Stratification and mixing (A 623; P12c), hydrostatic semi-implicit NHFLOW with density (A 620):
//
//  - the buoyancy frequency N^2 = -(g/rho0) d rho/dz at the layer centres, rho of the tracers
//    (nhflow_momentum_SI::eos), centred between the layers below and above (one-sided at the bed and the
//    surface), 0 in dry columns and with one layer
//  - the gradient Richardson number Ri = N^2/S^2, S^2 = (du/dz)^2 + (dv/dz)^2 the same way
//  - the stability functions of the algebraic eddy viscosity (A 547), Munk and Anderson (1948) by default:
//    nu = nu_n (1 + a_m Ri)^p_m, K = nu_n/sigma_t (1 + a_h Ri)^p_h (A 625 a_m p_m a_h p_h, default
//    10 -0.5 3.33 -1.5), Ri < 0: 1 (no enhancement)
//  - k-epsilon (A 560 1, 21): the buoyancy production G_b = -(nu_t/sigma_t) N^2 in the k and the epsilon
//    equations (nhflow_strain::Pk_strat, nhflow_kepsilon_func), epsilon with c1e c3e eps/k G_b, c3e = 1 for
//    G_b > 0 (unstable), c3e = A 624 for G_b < 0 (stable; default 0, Rodi 1987). Homogeneous steady state:
//    flux Richardson number -G_b/P = (c2e - c1e)/(c2e - c3e) = 0.25 for c3e = 0, Ri = sigma_t times that
// ---------------------------------------------------------------------------

void nhflow_n2(lexer*, fdm_nhf*, double *N2);
double nhflow_ri(lexer*, fdm_nhf*, const double *N2, int i, int j, int k);
double nhflow_fm(lexer*, double Ri);
double nhflow_fh(lexer*, double Ri);

#endif
