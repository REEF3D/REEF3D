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
--------------------------------------------------------------------
Parts of seastate_source.cpp are C++ translations of routines of
SWAN 41.51: FAC4WW and SWSNL1 (DIA quadruplets), FAC3WW and SWLTA
(LTA triads), SWIND5 (Yan wind input) and SWCAP with IWCAP 7 (van der
Westhuysen whitecapping), SSURF (Newton linearisation of the
Battjes-Janssen dissipation), SINTGRL (maximum energy with
Battjes-Janssen breaking) and SVEG (vegetation, one layer).
SWAN: Copyright (C) 1993-2024 Delft University of Technology, GNU GPL
version 3 or later (see seastate_source.cpp).
--------------------------------------------------------------------*/

#ifndef SEASTATE_SOURCE_H_
#define SEASTATE_SOURCE_H_

#include<vector>
#include <cstddef>
#include"seastate_st4.h"

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
                   (A 729 1: C_d of the met drag law W 207, met_drag.h, e.g.
                   capped in hurricane winds as the surge; A 729 2: the
                   law A 728 of its own, e.g. Zijlema for the waves)
  whitecapping     Komen et al. (1984) (SWCAP): C_ds 2.36e-5,
                   s_PM^2 3.02e-3, p 4, delta 1, mean sig_-10 and
                   k_WAM, spectral tail sig^-4                       -> D
  wind input and   van der Westhuysen et al. (2007), A 732 2 (SWAN GEN3
  whitecapping     WESTH): exponential growth of Yan (1987) (SWIND5),
  (A 732 2)        B = max(0, (0.04 s^2 + 0.00552 s + 0.000052) cos(theta-
                   theta_w) - 0.000302) sig, s = u_star/c            -> P
                   saturation-based whitecapping (SWCAP, IWCAP 7):
                   B_k = c_g k^3 E(sig), f_br = (1 + tanh(10 (sqrt(B_k/B_r)
                   - 1)))/2, p = 3 + tanh(25.76 (u_star k/sig - 0.1)),
                   f_br C_ds (B_k/B_r)^(p/2) (sqrt(g k)/sig)^(p/2-1)
                   sqrt(g k) + (1 - f_br) 3e-5 (st/s_PM)^4 sig_-10 k/k_WAM,
                   st = k_WAM sqrt(E_tot), C_ds 5e-5, B_r 1.75e-3 (A 757) -> D
  wind input,      Ardhuin et al. (2010) as the ECMWF model ecWAM, IPHYS 1
  swell damping,   (A 732 3, Phase 20, seastate_st4): Janssen growth rate with
  whitecapping     the sheltering of the friction velocity, u* and z0 from the
  (A 732 3)        stress balance with the wave stress of the spectrum (the
                   fixed point per cell), swell damping of Ardhuin et al.
                   (2009), saturation-based whitecapping B_r 9e-4, the
                   diagnostic f^-5 tail above max(2.5 f_m, 3 f_PM) and the
                   upper bound alpha_p 0.031 (see seastate_st4.h)    -> P, D
  quadruplets     DIA (Hasselmann et al. 1985) (FAC4WW, SWSNL1):
                   lambda 0.25, C 3e7, shallow-water scaling
                   1 + 5.5/x (1 - 0.833 x) exp(-1.25 x),
                   x = max(0.75 k_WAM d, 0.5)                       -> P, D
  breaking (2)     Roelvink (1993) for the wave groups of the surfbeat
                   mode (XBeach 'roelvink2'): Qb = 1 - exp(-(H/(gamma h))^n),
                   D/E = 2 alpha f_rep Qb H/h, H = sqrt(8 E)          -> D
  breaking (1)     Battjes and Janssen (1978) (SSURF, FRABRE):
                   alpha 1, gamma 0.73, mean frequency sig_01,
                   Newton linearisation as SWAN (SbrD)               -> P, D
  bottom friction  JONSWAP (SBOT): C_b/g^2 (sig/sinh(kd))^2,
                   C_b 0.038 m^2/s^3                                 -> D
  vegetation (1)   Dalrymple et al. (1984), Mendez and Losada (2004),
                   as SWAN (SVEG, IVEG 1, Suzuki et al. 2011), one layer:
                   sqrt(2/pi) g^2 (k/sig)^3 sqrt(E_tot) Cd bv Nv
                   (sinh^3 k ah + 3 sinh k ah)/(3 k cosh^3 kd),
                   k = k_WAM, sig = sig_01, ah = min(height, d),
                   zero for kd > 10                                  -> D
  vegetation (2)   per frequency, Jacobsen et al. (2019), as SWAN
                   (SVEG, IVEG 2): sqrt(2/pi)/g Cd bv Nv int_0^ah
                   S_u(sig, z) sqrt(mu_0(z)) dz, S_u = (sig cosh(kz)/
                   sinh(kd))^2 E, mu_0 = int S_u dsig, Simpson with 20
                   intervals                                         -> D
  triads           LTA of Eldeberky (1996), as the original SWAN LTA
                   (SWLTA): alpha_EB 0.05, sum frequencies below 2.5
                   sig_01, Madsen and Sorensen (1993) interaction
                   coefficient, biphase pi/2 (tanh(0.63/Ur) - 1),
                   active for Ursell numbers Ur >= 0.1               -> P, D
                   (cutfr, urcrit, urslim: A 736; the SWAN versions
                   40.x-41.31 used urcrit 0.2, urslim 0.01)
  maximum energy   with Battjes-Janssen breaking the total energy of a
                   cell is limited to (gamma d)^2/4 before the source
                   terms are evaluated (SWAN SINTGRL, A 737 1): cap

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
    int drag = 0;                        // drag of u_star: 0 Wu (1982), else the law of met_drag.h (A 729 1: W 207, A 729 2: A 728)
    double drag_cap = 2.5e-3, drag_floor = 1.5e-3;      // ... its cap and floor (W 213, W 214)

    bool komen = false;                  // whitecapping (and exponential wind input with wind)
    bool westh = false;                  // ... as van der Westhuysen et al. (2007): Yan wind input, saturation-based whitecapping (SWAN GEN3 WESTH)
    double cds2 = 5.0e-5, br = 1.75e-3;  // ... its coefficient and saturation threshold
    bool st4 = false;                    // ... as Ardhuin et al. (2010), ecWAM IPHYS 1 (A 732 3, Phase 20): wind input with sheltering,
                                         // swell damping, saturation-based whitecapping, u* from the wave stress, diagnostic tail
    seastate_st4_param st4p;
    bool dia = false;                    // quadruplets
    double limiter = 0.1;                // action density limiter gamma with komen, 0: off

    bool breaking = false;
    bool emax = false;                   // Battjes-Janssen: total energy limited to (gamma d)^2/4 (SWAN SINTGRL)
    int breaking_model = 1;              // 1 Battjes-Janssen, 2 Roelvink (1993) for wave groups (surfbeat)
    double alpha = 1.0, gamma = 0.73, nroel = 10.0;

    bool friction = false;
    double Cb = 0.038;

    bool triads = false;
    double alphaEB = 0.05, cutfr = 2.5, urcrit = 0.63, urslim = 0.1;

    int vegetation = 0;                  // 1 Dalrymple / Suzuki et al. (2011), 2 per frequency, Jacobsen et al. (2019)
    double vh = 0.0, vd = 0.0, vn = 0.0, vcd = 1.0;     // height [m], stem diameter [m], stems per m^2, drag coefficient

    bool any() const {return wind || komen || dia || breaking || friction || triads || vegetation>0;}
};

class seastate_source
{
public:
    seastate_source(const seastate_grid&, const seastate_source_param&);

    const seastate_source_param &param() const {return prm;}

    // wind field (A 730 2): wind of the cell before compute [m/s], direction the wind blows to [rad]
    void set_wind(double U10, double wdir) {prm.U10 = U10; prm.wdir = wdir;}

    // vegetation field (A 756 1): stems per m^2 of the cell before compute
    void set_vegetation(double nv) {prm.vn = nv;}

    // the parameters of another source object of the same grid (Phase 12: the copies kept for the threads)
    void sync_param(const seastate_source &o) {prm = o.prm;}

    // spectral sparsity (A 795): energy may appear in any bin of the cell (wind input, DIA)
    bool fills_spectrum() const {return (prm.wind && prm.U10>0.0) || prm.dia;}
    // LTA triads: the lowest frequency offset (negative) that feeds a bin, 0 without triads
    int triad_reach() const {return (prm.triads && nsig>1) ? tsm1 : 0;}
    // the triads acted in the last compute (Ursell number above urslim)
    bool triads_active() const {return prm.triads && ursell>=prm.urslim;}
    // row sums (over the directions per frequency) of the spectrum of the next cap()/compute(), known
    // to the caller (spectral sparsity: the sums over the active ranges)
    void set_rows(const double *r) const {for(int l=0; l<nsig; ++l) row[l] = r[l]; rows_ready = true;}

    // N of one cell (nbin), k and cg of the cell (nsig); P and D of size nbin (overwritten)
    // directions ma..mb and frequencies la..lb only (the quadrant of a sweep, the band of the spectral
    // sparsity; mb, lb < 0: all). The integral parameters use the whole spectrum; P and D outside are not
    // set. With DIA the transfers are computed for all bins.
    void compute(const float *N, double depth, const float *k, const float *cg, double *P, double *D, int ma=0, int mb=-1,
                 int la=0, int lb=-1);

    // maximum energy (prm.emax, Battjes-Janssen breaking): scales the spectrum of a cell down to the
    // total energy (gamma d)^2/4 (tail included, as SWAN SINTGRL); true if it was scaled
    bool cap(float *N, double depth) const;

    // single terms as dN/dt [N/s], for tests and diagnostics (S of size nbin, overwritten)
    void quadruplets(const float *N, double depth, const float *k, double *S, double *dSdN = nullptr);
    void triads(const float *N, double depth, const float *k, const float *cg, double *S);

    // maximum change of N per iteration of the last call (per frequency), nullptr: no limiter
    const double *limit() const {return (prm.komen && prm.limiter>0.0) ? dNmax.data() : nullptr;}

    // ST4 (A 732 3): keep u*, the growth rates and the cut-off of the last compute for the next computes of the
    // same cell (the widened windows of the spectral sparsity, A 795), true until freeze(false)
    void freeze(bool f) {frozen = f;}
    // ST4: the diagnostic tail after the solve of a cell, directions ma..mb (the cut-off of the last compute)
    void project(float *N, const float *k, const float *cg, int ma, int mb) const;
    // ST4: the first frequency of the diagnostic tail of the last compute, nsig without a tail
    int tail_from() const {return (prm.st4 && prm.wind && prm.U10>0.0) ? st4s.tail_from() : nsig;}
    // ST4: the upper bound of N per frequency (ecWAM FLMAX, Phillips alpha_p,max / pi per bin), nullptr: none
    const double *upper() const {return prm.st4 ? st4s.upper() : nullptr;}
    const seastate_st4 &st4() const {return st4s;}

    // integral parameters of the last call
    double Etot = 0.0, sigm01 = 0.0, sigm_10 = 0.0, km_wam = 0.0, Hs = 0.0, Qb = 0.0, ursell = 0.0;
    double ustar = 0.0;
    double brk_rate = 0.0;              // Roelvink breaking: dissipation rate D/E [1/s] of the last call (roller source)

    static double Qb_bj(double Hrms, double Hm);      // Battjes-Janssen fraction of breaking waves
    static double ustar_wu(double U10);               // friction velocity, drag of Wu (1982)
    double ustar_of(double U10) const;                // friction velocity of the drag prm.drag

private:
    void moments(const float *N, double depth, const float *k);
    void dia(const float *N, double depth);
    void lta(const float *N, double depth, const float *k, const float *cg);
    void vegetation(const float *N, double depth, const float *k, double *D);
    void split(const float *N, double *P, double *D);
    seastate_st4 st4s;
    bool frozen = false;
    int wa = 0, wb = 0;             // direction window of the current compute
    int fa = 0, fb = 0;             // frequency band of the current compute

    // DIA arrays with the directions continued periodically beyond 0 and ndir-1 (no index wrapping in
    // the loops): UE for m = -2 dpad .. ndir-1+2 dpad, SA1/SA2 and the derivatives for m = -dpad .. ndir-1+dpad
    double &ue(int l, int m) {return UE[size_t(l+uoff)*uw + size_t(m+2*dpad)];}
    double &sa1(int l, int m) {return SA1[sx(l,m)];}
    double &sa2(int l, int m) {return SA2[sx(l,m)];}
    size_t sx(int l, int m) const {return size_t(l+soff)*sw + size_t(m+dpad);}
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
    int dpad = 0, uw = 0, sw = 0;      // largest direction offset of the interactions, row lengths of UE and SA
    double awg[8], dal1, dal2, dal3;
    std::vector<double> UE, SA1, SA2, DA1C, DA1P, DA1M, DA2C, DA2P, DA2M, af11;

    // LTA: interpolation at sig/2 and 2 sig, E(f) and SA of one direction, coefficient per frequency
    int tsm, tsm1, tsp, tsp1;
    double twm, twm1, twp, twp1;
    std::vector<double> EL, SAL, TQ;
    std::vector<double> VE, VS, VQ;     // vegetation per frequency: E dsig, the integral over the height, S_u/E
    mutable std::vector<double> row;    // sums over the directions per frequency (rows)
    mutable bool rows_ready = false;
    double pw_smax = 0.0, pw_se0 = 0.0, pw_se1 = 0.0, pw_se2 = 0.0, pw_fachfr = 0.0;   // tail powers
    void rows(const float *N) const;
};

#endif
