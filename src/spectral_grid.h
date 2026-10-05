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

#ifndef SPECTRAL_GRID_H_
#define SPECTRAL_GRID_H_

#include<string>
#include<vector>

/*--------------------------------------------------------------------
REEF3D::Spectral - spectral (frequency-direction) grid

  frequencies: logarithmic, f_l = fmin * r^l, l = 0..nsig-1,
               r = (fmax/fmin)^(1/(nsig-1)); sig = 2 pi f [rad/s]
               bin widths dsig from the geometric midpoints between
               neighbouring frequencies (end bins mirrored)
  directions:  full circle, theta_m = m * dtheta, m = 0..ndir-1,
               dtheta = 2 pi / ndir; Cartesian convention: direction
               of propagation, counter-clockwise from the +x axis

  bin index of (l,m): l*ndir + m  (direction fastest), the layout of
  one spectrum in spectral_store

No lexer or MPI dependency, so the class is unit-testable stand-alone.
--------------------------------------------------------------------*/

class spectral_grid
{
public:
    spectral_grid(int nsig, double fmin, double fmax, int ndir);

    bool valid() const {return error.empty();}
    const std::string &message() const {return error;}

    int bin(int l, int m) const {return l*ndir + m;}

    int nsig, ndir, nbin;
    double fmin, fmax, ratio, dtheta;

    std::vector<double> f, sig, dsig;           // size nsig
    std::vector<double> theta, costh, sinth;    // size ndir

private:
    std::string error;
};

#endif
