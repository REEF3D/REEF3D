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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SCOUR_CORE_H_
#define SCOUR_CORE_H_

#include "increment.h"
#include "slice4.h"
#include <fstream>
#include <vector>

class lexer;
class slice;
class ghostcell;

using namespace std;

// P 187 1: empirical local scour and scour protection parameters from the near-bed velocity of FNPF
// (fnpf_scour) or NHFLOW (nhflow_scour), written with the bed output at the free-surface output
// frequency; P 198 x y: scour gauges.
//
// The solver fills ub, vb with the near-bed velocity of the undisturbed flow every time step (the
// structure itself is not resolved, P 188: diameter D), and for a computed current (NHFLOW) ucx, ucy
// with the depth-averaged velocity. The statistics are accumulated from P 188 t_s on (time weighted):
//   U_m   = sqrt(2) sigma_u, sigma_u^2 = var(u_b) + var(v_b)   (Sumer & Fredsoe 2001, irregular waves)
//   KC    = U_m T / D, T = P 188 T, or the wave input period (Tp of a spectrum, else T)
//   U_c   ambient current, P 189 (constant) or P 196 1 (file scour_current.dat: x y Ux Uy),
//         plus the time mean of the computed depth-averaged current (NHFLOW)
//   U_cw  = U_c/(U_c + U_m)
//   S/D   = 1.3 [1 - exp(-A (KC - B))], A = 0.03 + 0.75 U_cw^2.6, B = 6 exp(-4.7 U_cw)
//           (Sumer & Fredsoe 2002, live bed, tested for KC ~ 5 - 30), 0 for KC < B, 1.3 without waves
//   tau_c, tau_w, tau_m, tau_max: Soulsby (1997), rough bed z0 = d50/12, f_w = max(rough, smooth)
//   theta_max, theta_cr (Soulsby-Whitehouse) of the bed sediment (S 20 d50, S 22 density)
//   armour (scour protection, Soulsby procedure as in REEF_Tools; P 197 d rho_s alpha k_t gamma_E n_r):
//           U_m and U_c amplified by alpha k_t; k_s = n_r D50; f_w = 1.39 (A_m/(k_s/30))^-0.52,
//           A_m = U_m/omega; tau_w = 0.5 rho f_w U_m^2; current drag from Manning n = (2.5 D50)^(1/6)/20,
//           Chezy C = h^(1/6)/n, tau_c = rho g/C^2 U_c^2 (U_c depth-averaged); tau_m = tau_c [1 + 1.2
//           (tau_w/(tau_c+tau_w))^3.2]; tau_max = tau_m + tau_w (collinear, phi = 0); gamma_E tau_max
//           against tau_cr = 0.055 (rho_s - rho) g D50:
//           armour_SR = gamma_E tau_max/tau_cr of the stone d (> 1 not statically stable),
//           armour_d_req the D50 with gamma_E tau_max = tau_cr, armour_d_req_max its maximum over time
//           (accumulated over the bed outputs and, with P 199 dt > 0, every dt; at the gauges over the
//           gauge outputs; from 5 wave periods of statistics on)
// P 193 x y: scour protection gauges (REEF3D_<solver>_ScourProtection), U_m, U_c, the stresses of the
// stone d, armour_SR, the required D50 and its maximum over time
// P 199 dt Tw m: Tw > 0 exponentially weighted statistics (time scale Tw) instead of all from t_s on,
// for time-varying sea states; dt the output interval of the scour gauges P 198 x y; m the definition
// of U_m: 0 sqrt(2) sigma_u (irregular waves), 1 the maximum orbital velocity, from the extremes of
// u_b and v_b about their mean (regular waves, the crest velocity of nonlinear waves)
// (REEF3D_FNPF_Scour or REEF3D_NHFLOW_Scour, one file per gauge, the values of the cell that contains the gauge)

class scour_core : public increment
{
public:
    // solver: "FNPF" or "NHFLOW" (output folders); depth: the still water depth;
    // computed_current: ucx, ucy are filled by the solver
    scour_core(lexer *, ghostcell *, const char *solver, slice *depth, bool computed_current);
    virtual ~scour_core();

    // every time step, after the solver filled ub, vb (and ucx, ucy): the statistics
    void sample(lexer *, ghostcell *);

    // before a bed output: the scour parameters from the statistics
    void evaluate(lexer *, ghostcell *);

    // P 198: the scour gauges, every time step or every P 199 dt
    void gauges(lexer *, ghostcell *);  // and the scour protection gauges P 193
    void gauge_flush();

    // the output quantities: their names and fields
    static const int NQ = 15;
    static const char *name(int);
    slice *field(int);

    // bed velocity (instantaneous), and the derived fields
    slice4 ub, vb;
    slice4 Um, ubmax, KC, Uc, Ucw, SD, S;
    slice4 thmax, thcr, mob;
    slice4 armSR, armdreq, armdmax;

    // computed current of the solver (depth-averaged velocity), instantaneous
    slice4 ucx, ucy;

private:
    void current_file(lexer *, ghostcell *);
    double period(lexer *);
    void cell(lexer *, int, int, double, double *);
    void gauge_ini(lexer *);

    // Soulsby (1997): the maximum bed shear stress of waves and current
    double taumax(double Uc, double Um, double T, double h, double cosphi, double d) const;
    // Soulsby & Whitehouse (1997): the critical Shields parameter
    double thetacr(double d, double s) const;

    // scour protection (REEF_Tools): the stresses for the stone size d, and the required D50
    struct protection
    {
        double tw, tc, tm, tmax, tcr, sr;
    };
    protection protect(double Um, double Uc, double T, double h, double d) const;
    double protect_dreq(double Um, double Uc, double T, double h, double d0) const;
    void pgauge_ini(lexer *);
    void pgauges(lexer *, ghostcell *);
    void accumulate(lexer *);

    slice4 Ucx, Ucy;
    slice4 su, sv, suu, svv, suv;
    slice4 suc, svc;
    slice4 uhi, ulo, vhi, vlo;  // extremes of the bed velocity from t_s on (P 199 m 1)

    const char *solver;
    slice *depth;
    bool has_current;

    double wsum;
    double g, rho, nu;
    double p_aamp, p_kt, p_gE, p_nr, p_rhos;  // P 197 scour protection factors, stone density
    bool warned_T;

    // P 198 gauges
    int ngauge = 0;
    int gflush = 0;
    double nextgauge = 0.0;
    std::vector<int> gi, gj, gflag;
    std::vector<double> gval;
    std::vector<std::ofstream *> gout;
    std::vector<double> gmax;  // the maximum required D50 at the P 198 gauges

    // P 193 scour protection gauges
    static const int NP = 10;
    // the maximum over time collects from NACC wave periods of statistics on (the mean of the bed
    // velocity, and with it U_m, is not converged before)
    static constexpr double NACC = 5.0;
    int npg = 0;
    std::vector<int> pi_, pj_, pflag;
    std::vector<double> pval, pmax;
    std::vector<std::ofstream *> pout;
};

#endif
