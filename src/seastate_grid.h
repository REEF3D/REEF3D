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

#ifndef SEASTATE_GRID_H_
#define SEASTATE_GRID_H_

#include<string>
#include<vector>

/*--------------------------------------------------------------------
REEF3D::SEASTATE - spectral (frequency-direction) grid

  frequencies: logarithmic, f_l = fmin * r^l, l = 0..nsig-1,
               r = (fmax/fmin)^(1/(nsig-1)); sig = 2 pi f [rad/s]
               bin widths dsig from the geometric midpoints between
               neighbouring frequencies (end bins mirrored)
  directions:  full circle, theta_m = m * dtheta, m = 0..ndir-1,
               dtheta = 2 pi / ndir; Cartesian convention: direction
               of propagation, counter-clockwise from the +x axis

  bin index of (l,m): l*ndir + m  (direction fastest), the layout of
  one spectrum in seastate_store

  single-frequency grid (surfbeat, A 770 1): seastate_grid(frep, ndir),
               one frequency frep with dsig = 1, so that sig N(theta) is
               the frequency-integrated energy per direction E(theta)
               [m^2/rad] and the integrals of seastate_param apply
               unchanged (m0 = sum E dtheta)

  direction faces theta_m + dtheta/2 (costhf, sinthf) for the fluxes in
  theta; quad[m]: quadrant of direction m for the four sweeps of the
  implicit solver

  fine direction sector (A 715, sector()): the directions with centres in
  the sector th1 .. th2 (counter-clockwise) are divided into k bins each;
  the directions are then sorted by angle in [0, 2 pi), dth[m] is the
  width of bin m, its upper face theta_m + dth/2 (costhf, sinthf), quad[m]
  from the angle. wth[m] = dth[m]/dtheta, exactly 1 on the uniform grid,
  so that the integrals sum N wth and multiply by dtheta in both cases

No lexer or MPI dependency, so the class is unit-testable stand-alone.
--------------------------------------------------------------------*/

class seastate_grid
{
public:
    seastate_grid(int nsig, double fmin, double fmax, int ndir);
    seastate_grid(double frep, int ndir);            // single frequency (surfbeat)

    // fine direction sector (A 715): th1, th2 [rad], k bins per direction of the sector
    void sector(double th1, double th2, int k);

    // the bin that holds direction th [rad]
    int direction_bin(double th) const;

    bool valid() const {return error.empty();}
    const std::string &message() const {return error;}

    int bin(int l, int m) const {return l*ndir + m;}

    int nsig, ndir, nbin;
    double fmin, fmax, ratio, dtheta;

    std::vector<double> f, sig, dsig;           // size nsig
    std::vector<double> theta, costh, sinth;    // size ndir
    std::vector<double> costhf, sinthf;         // faces theta_m + dtheta/2, size ndir
    std::vector<int> quad;                      // quadrant 0..3 of direction m: theta in [q pi/2, (q+1) pi/2)
    std::vector<double> dth, wth;               // bin widths, dth/dtheta (1 on the uniform grid)
    bool uniform = true;

private:
    void directions();
    std::string error;
};

#endif
