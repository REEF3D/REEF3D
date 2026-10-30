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
seastate_st4.cpp is a C++ translation, changed for the SEASTATE
spectral grid and source term interface, of routines of the ECMWF
wave model ecWAM 1.5.14 (IPHYS = 1): SINPUT_ARD, SDISSIP_ARD,
TAUT_Z0, STRESSO, TAU_PHI_HF, FRCUTINDEX, FEMEANWS, FKMEAN (mean
frequency), CHNKMIN, INIT_X0TAUHF, INIT_SDISS_ARDH, TABU_SWELLFT,
KERKEI, KZEONE, IMPHFTAIL and the parameters of SETWAVPHYS.
ecWAM: Copyright 1996- ECMWF, licensed under the Apache License,
Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0); see
seastate_st4.cpp for the notice and the list of changes.
--------------------------------------------------------------------*/

#ifndef SEASTATE_ST4_H_
#define SEASTATE_ST4_H_

#include<vector>

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - wind input and dissipation of Ardhuin et al. (2010)
(WAVEWATCH III ST4) as in ecWAM IPHYS = 1, A 732 3 (Phase 20)

  wind input       Janssen (1991) growth rate with the sheltering of
                   Ardhuin et al. (2010): the friction velocity of a
                   frequency reduced by the wave stress of the lower
                   frequencies (TAUWSHELTER 0.25), beta_max 1.40,
                   z_alpha 0.008, rho_a/rho_w 1.225/1000      -> P = gamma N
  swell damping    Ardhuin et al. (2009): viscous and turbulent
                   boundary layer of the air over the swell, the
                   transition at the Reynolds number of the orbital
                   motion (SWELLF 0.66 ... SWELLF7 3.6e5)       -> D
  whitecapping     saturation-based (Ardhuin et al. 2010, TEST473):
                   B(k,theta) = k^3 c_g int E cos^2 over +-80 deg,
                   C_ds 2.2e-5 sig [0.3 (B_max/B_r - 1)^2 +
                   0.7 (B/B_r - 1)^2]_+, B_r 9e-4; no cumulative and
                   no wave-turbulence term (ecWAM SSDSC3 = SSDSC5 = 0) -> D

  friction velocity and roughness from the stress balance (TAUT_Z0):
  U10 = (ustar / kappa) ln(10/z0), z0 = alpha tau/sqrt(1 - tau_w/tau) / g +
  0.11 nu/ustar, the Charnock parameter capped for strong winds (CHNKMIN,
  alpha 0.0065 to 0.0001 around 33 m/s); the wave stress tau_w is the
  input of the prognostic spectrum up to the cut-off and of the f^-5
  tail above it (STRESSO, TAU_PHI_HF). ecWAM lags tau_w by a time step;
  here the fixed point ustar(tau_w(ustar)) of the spectrum of the cell is
  iterated from tau_w = 0 in every evaluation (no state, the same in any
  order of the cells), relative change of ustar below tol

  diagnostic tail (IMPHFTAIL): above the cut-off frequency
  f_c = max(2.5 max(f_m, f_m,ws), 3 g/(2 pi 28 ustar)) the spectrum is held
  at the saturation of f_c, E ~ 1/(k^3 c_g): after each solve of a cell
  its bins above the cut-off of the source terms of that solve are set
  to the tail (as ecWAM after each time step)

Units as seastate_source: N(sig,theta), P [N/s], D [1/s].
--------------------------------------------------------------------*/

struct seastate_st4_param
{
    double betamax = 1.40, zalp = 0.008;
    double alpha = 0.0065, alphamin = 0.0001, chnkmin_u = 33.0;
    bool capchnk = true;
    double tauwshelter = 0.25;
    bool swell = true;
    double swellf4 = 1.5e5, swellf7 = 3.6e5, z0rat = 0.04, z0tubmax = 0.0005;
    bool tail = true;
    double tailfactor = 2.5, tailfactor_pm = 3.0;
    double rhoa = 1.225, rhow = 1000.0, nu = 1.5e-5;
    double sdsbr = 9.0e-4, ssdsc2 = -2.2e-5, ssdsc4 = 1.0, ssdsc6 = 0.3;
    double tol = 1.0e-4;
    int maxit = 12;
    double alphapmax = 0.031;           // upper bound of the spectrum (ecWAM ALPHAPMAX), 0 off
};

class seastate_st4
{
public:
    void init(const seastate_grid&, const seastate_st4_param&, double gravity);
    const seastate_st4_param &param() const {return prm;}

    // the friction velocity and roughness for the spectrum N of a cell (rows: sums of N over the directions
    // per frequency), the input growth rates of all bins, the cut-off; freeze: those of the last call kept
    void wind(const float *N, const double *row, const float *k, const float *cg, double U10, double wdir, bool freeze);

    // no wind (a calm cell): no input and no swell damping, the whitecapping alone in add()
    void calm() {calmed = true; ust = 0.0; tauw = 0.0; lc = nsig-1;}

    // P, D of the input, swell damping and whitecapping for the frequencies la..lb and directions ma..mb
    // (after wind() of the same cell)
    void add(const float *N, const float *k, const float *cg, int la, int lb, int ma, int mb, double *P, double *D);

    // the diagnostic tail (IMPHFTAIL) after the solve of a cell: N above the cut-off of the last wind() set to the
    // saturation of the cut-off, directions ma..mb
    void project(float *N, const float *k, const float *cg, int ma, int mb) const;
    int tail_from() const {return (!prm.tail || calmed || !(ust>0.0)) ? nsig : lc+1;}

    // single terms of the last wind() for tests: growth rate gamma (1/s), swell damping (>= 0), whitecapping (>= 0)
    double rate(int b) const {return gam[b];}

    // the upper bound of N per frequency: F(f,theta) <= alpha_p,max/pi g^2 (2 pi)^-4 f^-5 (ecWAM FLMAX), nullptr: none
    const double *upper() const {return prm.alphapmax>0.0 ? nmax.data() : nullptr;}
    void dissipation(const float *N, const float *k, const float *cg, int l, double *Dl);

    // the state of the last wind()
    double ust = 0.0, z0 = 0.0, z0b = 0.0, chrnck = 0.0, tauw = 0.0, tauwdir = 0.0, fm = 0.0, fmws = 0.0;
    int lc = 0, iters = 0;
    bool calmed = true;

    static void kzeone(double x, double y, double &re0, double &im0, double &re1, double &im1);
    static void kerkei(double x, double &ker, double &kei);

private:
    void taut_z0(int iusfg, double U, double udir, double tw, double twdir, double &us, double &zz0, double &zz0b, double &chr) const;
    void input(const float *N, const float *k, const float *cg, double us, double zz0, double wdir, bool sneg);
    double tau_hf(int m, bool shelter, double us, double zz0, double wdir, const float *N, double ustp) const;
    void stress(const float *N, const float *k, double us, double zz0, double wdir, double &tw, double &twdir) const;
    void cutoff(const float *N, double us);
    double chnkmin(double u10) const;

    const seastate_grid *pg = nullptr;
    seastate_st4_param prm;
    int nsig = 0, ndir = 0, nbin = 0;
    double grav = 9.81, rnum = 0.0, raorw = 0.0, constn = 0.0;
    double x0tauhf = 0.0, wtauhf[19];
    std::vector<double> swellft;                // SWELLFT(1..200), index 0 unused
    int nsat = 0;                               // directions on each side of the saturation integral
    std::vector<int> satidx;                    // per direction 2 nsat + 1 indices, -1: no such direction (compact sector)
    std::vector<double> satw;
    std::vector<double> dfim, rhowg;            // DFIM [Hz rad] and RHOWG_DFIM of the SEASTATE grid
    std::vector<double> gam, dst;               // growth rate and swell damping (>= 0) per bin
    std::vector<double> ustl, udl, sxl, syl;    // ustar and its direction per frequency (sheltering), the stress of the input per frequency
    std::vector<double> wsl;                    // sum of F over the bins with positive input, per frequency
    std::vector<double> nmax;                   // upper bound of N per frequency
    std::vector<double> bth, wcl;               // saturation and whitecapping rate of one frequency
    double sw_temp2 = 0.0, sw_pturb = 0.0, sw_pvisc = 0.0;  // swell damping of the cell
};

#endif
