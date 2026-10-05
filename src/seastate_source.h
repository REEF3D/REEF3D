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

#ifndef SEASTATE_SOURCE_H_
#define SEASTATE_SOURCE_H_

#include<vector>

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - source terms of the wave action balance

For one cell, the source term of the action balance is returned split
as   dN/dt = P - D N   with P >= 0 (explicit) and D >= 0 (implicit), so
that the implicit transport matrix stays an M-matrix and N stays
non-negative without limiters (Patankar splitting; positive parts of
the nonlinear transfers go into P, negative parts into D). The
quadruplets are linearised around the latest spectrum with their
diagonal derivative L = dS/dN (as SWAN, SWSNL1): where L < 0,
S(N_new) ~ S(N) + L (N_new - N) is split as P += -L N, D += -L, which
keeps P, D >= 0 and damps the stiff high-frequency transfer (a damped
Newton step instead of a fixed-point iteration).

Action density limiter (SWAN, Ris 1997), with the deep-water physics
only (as SWAN GEN3 and WAM): the change of N per iteration is limited
to gamma alpha_PM/(2 sig k^3 c_g), gamma = A 735 (0.1), alpha_PM =
0.0081. It is active only while the wind-sea spectrum is far from its
balance (large time steps, early stationary iterations); breaking,
friction and triads are not limited.

The formulations and default coefficients follow SWAN 41.51 (routines
in brackets):

  wind input       linear growth (Cavaleri and Malanotte-Rizzoli 1981)
                   with the filter of Tolman (1992) (SWIND0), A 734;
                   exponential growth of Komen et al. (1984) (SWIND3),
                   B = max(0, 0.25 rho_a/rho_w (28 u_star cos(theta-theta_w)/c - 1)) sig,
                   u_star = sqrt(C_d) U10 with the drag of Wu (1982)     -> P
  whitecapping     Komen et al. (1984) (SWCAP): C_ds 2.36e-5,
                   s_PM^2 3.02e-3, p 4, delta 1, mean sig_-10 and
                   k_WAM, spectral tail sig^-4                       -> D
  quadruplets      DIA (Hasselmann et al. 1985) (FAC4WW, SWSNL1):
                   lambda 0.25, C 3e7, shallow-water scaling
                   1 + 5.5/x (1 - 0.833 x) exp(-1.25 x),
                   x = max(0.75 k_WAM d, 0.5)                       -> P, D
  breaking         Battjes and Janssen (1978) (SSURF, FRABRE):
                   alpha 1, gamma 0.73, mean frequency sig_01,
                   Newton linearisation as SWAN (SbrD)               -> P, D
  bottom friction  JONSWAP (SBOT): C_b/g^2 (sig/sinh(kd))^2,
                   C_b 0.038 m^2/s^3                                 -> D
  triads           LTA of Eldeberky (1996), as the original SWAN LTA
                   (SWLTA): alpha_EB 0.05, sum frequencies below 2.5
                   sig_01, Madsen and Sorensen (1993) interaction
                   coefficient, biphase pi/2 (tanh(0.63/Ur) - 1),
                   active for Ursell numbers Ur >= 0.1               -> P, D

Integral parameters (SINTGRL): E_tot, sig_01, sig_-10, k_WAM with a
sig^-4 tail above the frequency grid, Hs = 4 sqrt(E_tot),
Ur = g Hs / (2 sqrt(2) sig_01^2 d^2).

Units: N(sig,theta) [m^2 s^2/rad^2], E = sig N, P [N/s], D [1/s].
No lexer or MPI dependency (unit test Regression/unit/seastate_test.cpp).
--------------------------------------------------------------------*/

struct seastate_source_param
{
    bool wind = false;
    double U10 = 0.0, wdir = 0.0;        // [m/s], direction the wind blows to [rad]
    double Alin = 1.5e-3;                // linear growth coefficient, 0: off

    bool komen = false;                  // whitecapping (and exponential wind input with wind)
    bool dia = false;                    // quadruplets
    double limiter = 0.1;                // action density limiter gamma with komen, 0: off

    bool breaking = false;
    double alpha = 1.0, gamma = 0.73;

    bool friction = false;
    double Cb = 0.038;

    bool triads = false;
    double alphaEB = 0.05, cutfr = 2.5, urcrit = 0.63, urslim = 0.1;

    bool any() const {return wind || komen || dia || breaking || friction || triads;}
};

class seastate_source
{
public:
    seastate_source(const seastate_grid&, const seastate_source_param&);

    const seastate_source_param &param() const {return prm;}

    // N of one cell (nbin), k and cg of the cell (nsig); P and D of size nbin (overwritten)
    void compute(const float *N, double depth, const float *k, const float *cg, double *P, double *D);

    // single terms as dN/dt [N/s], for tests and diagnostics (S of size nbin, overwritten)
    void quadruplets(const float *N, double depth, const float *k, double *S, double *dSdN = nullptr);
    void triads(const float *N, double depth, const float *k, const float *cg, double *S);

    // maximum change of N per iteration of the last call (per frequency), nullptr: no limiter
    const double *limit() const {return (prm.komen && prm.limiter>0.0) ? dNmax.data() : nullptr;}

    // integral parameters of the last call
    double Etot = 0.0, sigm01 = 0.0, sigm_10 = 0.0, km_wam = 0.0, Hs = 0.0, Qb = 0.0, ursell = 0.0;
    double ustar = 0.0;

    static double Qb_bj(double Hrms, double Hm);      // Battjes-Janssen fraction of breaking waves
    static double ustar_wu(double U10);               // friction velocity, drag of Wu (1982)

private:
    void moments(const float *N, double depth, const float *k);
    void dia(const float *N, double depth);
    void lta(const float *N, double depth, const float *k, const float *cg);
    void split(const float *N, double *P, double *D);

    double &ue(int l, int m) {return UE[size_t(l+uoff)*ndir + m];}
    double &sa1(int l, int m) {return SA1[size_t(l+soff)*ndir + m];}
    double &sa2(int l, int m) {return SA2[size_t(l+soff)*ndir + m];}
    size_t sx(int l, int m) const {return size_t(l+soff)*ndir + m;}
    int dw(int m) const {return ((m%ndir)+ndir)%ndir;}

    const seastate_grid &g;
    seastate_source_param prm;
    int nsig, ndir, nbin;
    double lnr;

    std::vector<double> S;          // nonlinear transfers dN/dt, split into P and D
    std::vector<double> L;          // diagonal derivative dS/dN of the quadruplets
    std::vector<double> dNmax;      // action density limiter per frequency

    // DIA (FAC4WW): interpolation indices and weights, E(f,theta) incl. the tail, SA1/SA2
    int isp, isp1, ism, ism1, idp, idp1, idm, idm1, uoff, ulen, soff, slen;
    double awg[8], dal1, dal2, dal3;
    std::vector<double> UE, SA1, SA2, DA1C, DA1P, DA1M, DA2C, DA2P, DA2M, af11;

    // LTA: interpolation at sig/2 and 2 sig, E(f) and SA of one direction, coefficient per frequency
    int tsm, tsm1, tsp, tsp1;
    double twm, twm1, twp, twp1;
    std::vector<double> EL, SAL, TQ;
};

#endif
