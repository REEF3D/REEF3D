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

#ifndef SHIP_MODELS_H_
#define SHIP_MODELS_H_

#include<vector>

//  Semi-empirical ship load models (solver independent, no lexer). Ship frame: x forward,
//  y to port, z up, origin at the centre of gravity; u, v, w and p, q, r are the velocities
//  in this frame.
//
//  friction:   ITTC-1957 correlation line with form factor (ITTC 7.5-02-03-01.4),
//              X = -1/2 rho S (1+k) C_F |u| u,  C_F = 0.075/(log10 Re - 2)^2,  Re = |u| L / nu
//  crossflow:  strip theory cross-flow drag (Hooft 1994),
//              Y = -1/2 rho C_D sum T(x) |v + x r| (v + x r) dx,   N = the same with the lever x
//  roll:       linear and quadratic roll damping  K = -B44 p - B44q |p| p
//  propeller:  T = rho n^2 D^4 KT(J), Q = rho n^2 D^5 KQ(J), J = Va/(n D)
//  rudder:     MMG standard method (Yasukawa & Yoshimura 2015)

class ship_models
{
public:
    
    // ITTC-1957 friction coefficient; Re below Re_min (1e5, below the turbulent range of the
    // correlation line) is taken as Re_min
    static double cf_ittc57(double Re);
    static const double Re_min;
    
    // frictional resistance including the form factor; returns X, Re and C_F
    static double friction(double rho, double nu, double S, double L, double k, double u, double &Re, double &CF);
    
    // cross-flow drag of the strips (xs: strip centres relative to the CoG, dx: widths, T: drafts)
    static void crossflow(double rho, double Cd, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
                          double v, double r, double &Y, double &N);
    
    static double roll_damping(double B44, double B44q, double p);
    
    // propeller open-water characteristics KT(J) = kt0 + kt1 J + kt2 J^2, KQ(J) likewise;
    // n [rev/s], D [m], Va axial inflow speed: returns J, KT, KQ, thrust T and torque Q
    static void propeller(double rho, double n, double D, const double *kt, const double *kq, double Va,
                          double &J, double &KT, double &KQ, double &T, double &Q);
    
    // MMG standard rudder (Yasukawa & Yoshimura 2015), converted to the ship frame here (y to
    // port, z up). delta > 0: rudder to starboard (trailing edge to starboard), the ship turns to
    // starboard. Positions relative to the CoG.
    struct rudder_param
    {
        double AR;      // rudder area [m^2]
        double Lambda;  // aspect ratio
        double xR, zR;  // rudder position (xR < 0 aft) [m]
        double tR;      // steering resistance deduction
        double aH;      // rudder force increase factor
        double xH;      // position of the additional lateral force [m]
        double eps;     // wake ratio of propeller and rudder (1-wR)/(1-wP)
        double kappa;   // propeller slipstream factor
        double lR;      // effective longitudinal position of the rudder for the flow straightening [m]
        double gammaR;  // flow straightening coefficient
    };
    
    // u, v, r: ship velocities; delta [rad]; propeller (D, n, KT, uP = (1-wP) u; D = 0: none):
    // returns X, Y, N, K and the effective rudder angle alpha_R, inflow speed U_R, normal force F_N
    static void rudder_mmg(double rho, const rudder_param&, double u, double v, double r, double delta,
                           double D, double n, double KT, double uP,
                           double &X, double &Y, double &N, double &K, double &alphaR, double &UR, double &FN);
};

#endif
