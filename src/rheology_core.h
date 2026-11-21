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

#ifndef RHEOLOGY_CORE_H_
#define RHEOLOGY_CORE_H_

// Rheology core: the constitutive physics of non-Newtonian, granular and debris materials, without
// any dependence on the hydrodynamic solvers (no lexer, fdm, fdm2D, fdm_nhf). The solvers connect
// through thin adapters (SFLOW: sflow_rheology_core; CFD and NHFLOW: later stages).
//
// 1. local (3D) laws: effective dynamic viscosity mu_eff(gamma, p') for the momentum equations
//    - Herschel-Bulkley with Papanastasiou regularisation, yield stress constant or pressure
//      dependent (Herschel-Bulkley-Coulomb / Drucker-Prager)
//    - mu(I) (GDR MiDi 2004, Jop et al. 2006) with the partial regularisation of Barker & Gray
//      (2017) at low I (well-posed for I < I_N, I_N from the criterion of Barker et al. 2015)
//    - Coulomb-viscoplastic (Domnik et al. 2013)
// 2. depth-integrated closures: basal shear stress tau_b(h, U, sigma_n') for depth-averaged layers
//    - Coulomb, Voellmy, Pouliquen-Forterre mu(h,Fr), depth-averaged mu(I) (Bagnold profile),
//      Herschel-Bulkley exact plug-flow law (Bingham for n = 1), Herschel-Bulkley-Coulomb
//    - backward-Euler friction step with the yield criterion: the layer stops when the driving
//      force is below the static threshold (no sign chatter, rest on slopes below the friction angle)
// 3. basal pore pressure: effective normal stress sigma_n' = (1 - lambda) sigma_n, lambda relaxing
//    to its hydrostatic value with the diffusion time scale 4 h^2/(pi^2 D) (Iverson & Denlinger 2001)
//
// Units: SI, stresses in Pa, dynamic viscosities in Pa s, angles in radians inside the core.

namespace rheo
{

constexpr double PI_ = 3.14159265358979323846;

inline double deg2rad(double a) { return a*PI_/180.0; }

// ---------------------------------------------------------------------------------------------
// 1. local laws
// ---------------------------------------------------------------------------------------------

enum law3D
{
    LAW_NEWTONIAN    = 0,
    LAW_HB           = 1,   // Herschel-Bulkley, constant yield stress
    LAW_HBC          = 2,   // Herschel-Bulkley-Coulomb: tau_y = c + tan(phi) p'
    LAW_MUI          = 3,   // mu(I), partially regularised
    LAW_COULOMB_VP   = 4    // Coulomb-viscoplastic
};

struct law3D_params
{
    int law = LAW_NEWTONIAN;
    double mu0 = 1.0e-3;          // Newtonian viscosity [Pa s]
    double tauy = 0.0;            // yield stress / cohesion [Pa]
    double K = 1.0e-3;            // consistency [Pa s^n]
    double n = 1.0;               // exponent
    double tanphi = 0.0;          // tan of the friction angle (HBC, Coulomb-viscoplastic)
    double m = 1000.0;            // Papanastasiou parameter [s]
    double mu_max = 1.0e6;        // cap of the effective viscosity [Pa s]
    double mu_min = 0.0;          // floor of the effective viscosity [Pa s]
    // mu(I)
    double mus = 0.38;            // tan(delta_s)
    double mu2 = 0.64;            // tan(delta_2)
    double I0 = 0.279;
    double d = 5.0e-4;            // grain diameter [m]
    double rho_s = 2500.0;        // grain density [kg/m^3]
    double alpha = 1.9;           // partial regularisation parameter (alpha < 2: well-posed branch)
    double I_N = 0.0;             // switch to the regularised branch (set by mu_I_setup)
    double A = 0.0;               // continuity constant of the regularised branch (set by mu_I_setup)
};

// inertial number I = gamma d / sqrt(p'/rho_s)
double inertial_number(double gamma, double p_eff, double d, double rho_s);

// unregularised mu(I)
double mu_I_plain(const law3D_params &q, double I);

// partially regularised mu(I) (Barker & Gray 2017): sqrt(alpha/ln(A/I)) below I_N, mu(I) above
double mu_I(const law3D_params &q, double I);

// Barker et al. (2015): C = 4 nu^2 - 4 nu + mu^2 (1 - nu/2)^2 with nu = I mu'(I)/mu(I);
// C <= 0 well posed, C > 0 ill posed
double well_posedness(double mu, double nu);

// I_N (lower neutral point of the plain law, C = 0) and A (continuity at I_N); returns the upper
// neutral point I_2N (the plain law is ill posed again above it), or 0 when the plain law is ill
// posed everywhere
double mu_I_setup(law3D_params &q);

// yield stress [Pa] from the effective pressure (no tension: p' < 0 counts as 0)
double yield_stress(const law3D_params &q, double p_eff);

// effective dynamic viscosity [Pa s] from the shear rate sqrt(2 D:D) and the effective pressure
double viscosity(const law3D_params &q, double gamma, double p_eff);

// ---------------------------------------------------------------------------------------------
// 2. depth-integrated closures
// ---------------------------------------------------------------------------------------------

enum basal_law
{
    BASAL_OFF     = 0,
    BASAL_COULOMB = 1,
    BASAL_VOELLMY = 2,
    BASAL_PF      = 3,      // Pouliquen & Forterre (2002) mu(h,Fr)
    BASAL_MUI     = 4,      // mu(I) with the depth-averaged inertial number of the Bagnold profile
    BASAL_HB      = 5,      // Herschel-Bulkley plug flow (Bingham for n = 1)
    BASAL_HBC     = 6       // Herschel-Bulkley plug flow with tau_y = c + tan(delta) sigma_n'
};

struct basal_params
{
    int law = BASAL_OFF;
    double tandelta = 0.577;      // Coulomb, Voellmy, HBC
    double xi = 1000.0;           // Voellmy [m/s^2]
    // Pouliquen-Forterre
    double tz1 = 0.384, tz2 = 0.594, tz3 = 0.408;     // tan(zeta_1,2,3)
    double beta = 0.136, L = 6.5e-4;
    double kappa = 1.0e-3;        // exponent of the intermediate regime
    // mu(I)
    double mus = 0.38, mu2 = 0.64, I0 = 0.279, d = 5.0e-4, phis = 0.6;
    // Herschel-Bulkley
    double tauy = 0.0, K = 1.0, n = 1.0;
};

struct basal_state
{
    double h = 0.0;               // layer depth [m]
    double rho = 1000.0;          // bulk density [kg/m^3]
    double gn = 9.81;             // gravity normal to the bed [m/s^2]
    double g = 9.81;              // magnitude of gravity (Voellmy turbulent friction) [m/s^2]
    double sigma_n = 0.0;         // basal normal stress rho gn h [Pa]
    double lambda = 0.0;          // pore pressure ratio: sigma_n' = (1 - lambda) sigma_n
};

// basal shear stress magnitude [Pa] at the depth-averaged speed U >= 0
double basal_stress(const basal_params &q, const basal_state &s, double U);

// static threshold tau_b(U -> 0+) [Pa]
double basal_stress_static(const basal_params &q, const basal_state &s);

// Herschel-Bulkley plug flow: depth-averaged speed from the basal stress, and its inverse
double hb_mean_speed(double tau_b, double tauy, double K, double n, double h);
double hb_basal_stress(double U, double tauy, double K, double n, double h);

// backward-Euler friction step: the speed s >= 0 with s + dt tau_b(s)/(rho h) = Ustar; zero when
// Ustar <= dt tau_static/(rho h) (the layer stops or stays at rest)
double implicit_speed(const basal_params &q, const basal_state &s, double Ustar, double dt);

// ---------------------------------------------------------------------------------------------
// 3. pore pressure
// ---------------------------------------------------------------------------------------------

// relaxation of the pore pressure ratio over dt: lambda_h + (lambda - lambda_h) exp(-pi^2 D dt/(4 h^2))
double pore_relax(double lambda, double lambda_h, double D, double h, double dt);

}

#endif
