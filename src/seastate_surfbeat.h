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

#ifndef SEASTATE_SURFBEAT_H_
#define SEASTATE_SURFBEAT_H_

#include<vector>

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE surfbeat - wave-group boundary conditions (A 770 1)

From a frequency-direction spectrum (the SEASTATE boundary spectrum on
a multi-frequency grid) a random-phase wave train is built and turned
into the time series that drive the wave-group (surfbeat) model and
the long waves of the host model at the offshore boundary:

components   f_m = m df, df = 1/T_rec (A 772), over the range where the
             frequency spectrum exceeds 1e-3 of its peak; amplitude
             a_m = sqrt(2 S(f_m) df); direction drawn at random from
             the directional distribution of the spectrum at f_m
             (single summation, short-crested groups); random phase.
             Seed A 773: identical on all ranks and runs.

envelope     Z(t,y) = sum_m a_m exp(i(2 pi f_m t - k_m sin(theta_m) y + phi_m))
             E(t,y) = |Z|^2 / 2 [m^2]   (mean over T_rec = m0)
             E(theta,t,y) = E(t,y) Dbar(theta), Dbar the energy-weighted
             mean directional distribution (sum Dbar dtheta = 1)

bound wave   second-order difference interactions (Herbers et al. 1994,
             eq. A5, for surface elevation as XBeach, Van Dongeren et al.
             2003): for every pair m < n with df_mn = f_n - f_m < f_m,
               zeta_b = sum D_mn a_m a_n cos(dpsi_mn),
               dpsi = (w_n - w_m) t - (k_n - k_m).x + phi_n - phi_m,
             flux q_b = zeta_b c3 (cos theta3, sin theta3), c3 = dw/|dk|
             (limited to 0.8 of the free long-wave celerity, XBeach nmax),
             theta3 the direction of dk. In shallow water D tends to
             -g (2n - 1/2)/(g h - c_g^2) (Longuet-Higgins and Stewart 1962).

Time series are precomputed on dt_bc = T_rep/10 for every boundary row
y and repeat with the period T_rec. No lexer or MPI dependency
(unit test Regression/unit/seastate_test.cpp).
--------------------------------------------------------------------*/

class seastate_surfbeat
{
public:
    // gin: multi-frequency grid of the input spectrum Nin (action density)
    seastate_surfbeat(const seastate_grid &gin, const float *Nin, double trec, int seed);

    // Herbers (1994) difference-interaction coefficient for surface elevation (m^-1)
    static double herbers(double f1, double th1, double k1, double f2, double th2, double k2, double h, double &c3, double &th3);

    // time series for the boundary rows at y (the boundary is the line x = 0 of the components)
    void series(const std::vector<double> &y, const std::vector<double> &depth, bool bound);

    // values of row r at time t (periodic, linear in time)
    void at(int r, double t, double &E, double &zeta, double &qx, double &qy) const;

    int K = 0;                                      // number of components
    double df = 0.0, trep = 0.0, frep = 0.0, m0 = 0.0, dtbc = 0.0, trec = 0.0;
    std::vector<double> fk, ak, thk, phk;           // components
    std::vector<double> Dbar;                       // mean directional distribution on the directions of gin

private:
    const seastate_grid &g;
    int nt = 0;
    std::vector<std::vector<double>> Et, Zt, Qx, Qy;
};

#endif
