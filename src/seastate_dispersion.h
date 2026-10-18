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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SEASTATE_DISPERSION_H_
#define SEASTATE_DISPERSION_H_

/*--------------------------------------------------------------------
REEF3D::SEASTATE - linear dispersion relation sig^2 = g k tanh(k d)

  seastate_wavenumber(sig,d)  k: Newton iteration from the explicit
                              approximation of Guo (2002), relative
                              accuracy 1e-12; deep water (k d > 30 for
                              the deep-water k) k = sig^2/g
  seastate_cg(sig,k,d)        group velocity n sig/k,
                              n = 1/2 (1 + 2kd/sinh(2kd)), kd capped at 30
  seastate_refraction(sig,k,d) sig/sinh(2kd), the depth-refraction and
                              frequency-shift coefficient (0 for kd > 30)

No lexer dependency (unit test Regression/unit/seastate_test.cpp).
--------------------------------------------------------------------*/

const double seastate_gravity = 9.81;

double seastate_wavenumber(double sig, double d);
double seastate_cg(double sig, double k, double d);
double seastate_refraction(double sig, double k, double d);

#endif
