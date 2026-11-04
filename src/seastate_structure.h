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

#ifndef SEASTATE_STRUCTURE_H_
#define SEASTATE_STRUCTURE_H_

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - transmission formulas of the structures (A 725,
seastate_obstacle); no lexer or MPI dependency (unit test
Regression/unit/seastate_test.cpp)

  seastate_dangremond  d'Angremond, Van der Meer and De Jong (1996), as
                       SWAN DAM DANGREMOND: the wave-height transmission
                       Kt from the freeboard Rc, Hs, the peak period Tp,
                       the seaward slope [deg] and the crest width B
  seastate_porous      porous slab through the water column per
                       frequency (Madsen 1974, Sollitt and Cross 1972,
                       van Gent 1995 resistance, see seastate_obstacle.h):
                       Kt^2 and Kr^2 for the incident energy E(sig)
                       [m^2 s/rad] per frequency of grid g, depth h,
                       width B, porosity n, stone diameter D50; returns
                       the rms discharge velocity in the structure
--------------------------------------------------------------------*/

double seastate_dangremond(double Rc, double Hs, double Tp, double slope, double B);

double seastate_porous(const seastate_grid &g, const double *E, double h, double B, double n, double D50,
                       float *kt2, float *kr2);

#endif
