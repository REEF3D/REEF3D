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

#ifndef SEASTATE_PARAM_H_
#define SEASTATE_PARAM_H_

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - integrated wave parameters of one spectrum

  energy density  E(sig,theta) = sig * N(sig,theta)    [m^2 s / rad^2]
  moments         m_n = sum_l sum_m sig^n E dsig dtheta
  Hs     = 4 sqrt(m0)                 significant wave height (Hm0)
  Tm01   = 2 pi m0/m1                 mean period
  Tm10   = 2 pi m_-1/m0               energy period Te (Tm-1,0)
  Tp     = 2 pi / sig_peak            peak period of E(sig), discrete
  dir    = atan2(b,a) [deg, 0..360)   mean direction of propagation
                                      (Cartesian, ccw from +x)
  spread = sqrt(2(1 - sqrt(a^2+b^2)/m0)) [deg]  directional spread (Kuik 1988)
  with a = sum E cos(theta), b = sum E sin(theta)
--------------------------------------------------------------------*/

struct seastate_param
{
    double m0 = 0.0, m1 = 0.0, mm1 = 0.0;
    double Hs = 0.0, Tm01 = 0.0, Tm10 = 0.0, Tp = 0.0, dir = 0.0, spread = 0.0;
    int lpeak = -1;

    void compute(const seastate_grid &g, const float *N);
};

#endif
