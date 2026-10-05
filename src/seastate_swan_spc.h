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

#ifndef SEASTATE_SWAN_SPC_H_
#define SEASTATE_SWAN_SPC_H_

#include<string>
#include<vector>

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - SWAN standard spectral file (.spc), 2D spectra

  read     first location and first time of the file: frequencies
           (AFREQ or RFREQ, Hz), directions (CDIR Cartesian or NDIR
           nautical, converted to Cartesian direction of propagation:
           theta = 270 - nautical, as SWAN with the default north
           direction), QUANT VaDens [m2/Hz/degr] or EnDens [J/m2/Hz/degr]
           (divided by rho g = 1025 * 9.81), FACTOR/ZERO/NODATA blocks
  to_grid  bilinear interpolation of E(f,theta) onto the spectral grid
           (linear in f, zero outside the file's frequency range,
           periodic linear in theta) and conversion to action density
           N(sig,theta) = E(f,theta[deg]) * 180/pi / (2 pi) / sig

No lexer dependency (unit test Regression/unit/seastate_test.cpp).
--------------------------------------------------------------------*/

struct seastate_swan_spc
{
    std::vector<double> f;          // Hz, ascending
    std::vector<double> dir;        // Cartesian direction of propagation [deg], 0..360
    std::vector<double> E;          // variance density [m2/Hz/degr], E[nf*ndir]: f major
    double x = 0.0, y = 0.0;        // location

    bool read(const std::string &file, std::string &error);
    void to_grid(const seastate_grid &g, std::vector<float> &N) const;
};

#endif
