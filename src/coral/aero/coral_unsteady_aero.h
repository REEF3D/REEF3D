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
// REEF3D::CORAL: unsteady aerodynamics of a blade element: dynamic inflow and dynamic stall, as AeroDyn's DBEMT
// (Øye) and UA (UAMod 4: Beddoes-Leishman in the state-space form of Hansen, Gaunaa and Madsen 2004; UAMod 6: Øye).
//
// Dynamic inflow (Øye), for the induced velocity W = (W_n, W_t) of a station (W_n = a Vx along the rotor normal,
// reducing the inflow; W_t = a' Vy in the plane), from its quasi-steady value W_qs (blade element momentum):
//   W_int + tau1 dW_int/dt = W_qs + k tau1 dW_qs/dt,    W + tau2 dW/dt = W_int
//   tau1 = 1.1/(1 - 1.3 min(a_mean, 0.5)) R/U0,          tau2 = (0.39 - 0.26 (r/R)^2) tau1,   k = 0.6
// with a_mean the axial induction averaged over the rotor and U0 the wind speed at the hub.
//
// Dynamic stall, the inputs of a step: the angle of attack at the aerodynamic centre alpha_ac, the relative speed U
// and the section's nose-up angular velocity omega (pitch and torsion); T_u = c/(2U) in [0.001, 50] s,
// alpha_34 = alpha_ac + T_u omega (the 3/4-chord point). From the static polar (stall_polar):
//   Cl_alpha the largest slope of cl between -20 and 20 deg (or C_lalpha of the file), alpha0 the zero lift angle
//   (or alpha0 of the file), Cl_inv = Cl_alpha (alpha - alpha0);
//   f_st = min((2 sqrt(cl/Cl_inv) - 1)^2, 1) from alpha0 outwards, 0 beyond the first angle where 2 sqrt(cl/Cl_inv) <= 1
//   on each side, linear between the table's angles; Cl_fs = (cl - Cl_inv f_st)/(1 - f_st) (the fully separated lift)
//   from the interpolated f_st and the static cl at the angle, cl/2 where f_st = 1 (within 1e-3: the table's value).
// Øye (model 1): df/dt = (f_st(alpha_34) - f)/T_f, T_f = T_f0 T_u (T_f0 = 6);
//   cl = f Cl_inv(alpha_34) + (1 - f) Cl_fs(alpha_34), cd and cm static at alpha_ac.
// Beddoes-Leishman HGM (model 2), four states:
//   dx1/dt = b1/T_u (A1 alpha_34 - x1),  dx2/dt = b2/T_u (A2 alpha_34 - x2),  alpha_E = (1 - A1 - A2) alpha_34 + x1 + x2
//   Cl_p = Cl_alpha (alpha_E - alpha0) + pi T_u omega,  dx3/dt = (Cl_p - x3)/T_p,  alpha_F = x3/Cl_alpha + alpha0
//   dx4/dt = (f_st(alpha_F) - x4)/T_f,  T_p = T_p0 T_u, T_f = T_f0 T_u
//   Cl_circ = x4 Cl_alpha (alpha_E - alpha0) + (1 - x4) Cl_fs(alpha_E),   cl = Cl_circ + pi T_u omega
//   cd = cd(alpha_E) + (alpha_34 - alpha_E) Cl_circ + (cd(alpha_E) - cd(alpha0)) dCd,
//        dCd = (sqrt(f_st(alpha_E)) - sqrt(x4))/2 - (f_st(alpha_E) - x4)/4
//   cm = cm(alpha_E) - pi/2 T_u omega
//   A1 0.165, A2 0.335, b1 0.0455, b2 0.3 (Jones' approximation of Wagner's function, as Hansen, Gaunaa and Madsen;
//   its attached-flow response follows Theodorsen's function; AeroDyn's DEFAULT is 0.3, 0.7, 0.14, 0.53), T_f0 3,
//   T_p0 1.7 (AeroDyn's defaults), or the values of the airfoil file.
//
// Time integration: every state is a first-order lag dy/dt = (u - y)/tau. A step from the committed state at t0 to
// t1 = t0 + dt takes u linear in time between its values at t0 and t1 and integrates exactly (lag_step):
//   y1 = u1 + (y0 - u0) E - (u1 - u0) tau/dt (1 - E),  E = exp(-dt/tau)
// stable for any dt/tau, exact for a linear input; tau is held over the step. The values at t0 are kept with the
// states, so a step only needs the committed state and the inputs at t1: the trial states of the implicit time
// integrator give loads consistent with the step, and the converged step commits them.

#ifndef CORAL_UNSTEADY_AERO_H_
#define CORAL_UNSTEADY_AERO_H_

#include"coral/aero/coral_airfoil.h"
#include<vector>

namespace coral
{
    // one step of dy/dt = (u - y)/tau with u linear from u0 to u1 over dt (exact); dt <= 0: y0
    double lag_step(double y0, double u0, double u1, double tau, double dt);

    struct dynamic_inflow_settings
    {
        bool on = false;
        double k = 0.6;
        double tau1 = 0.0;          // > 0: fixed tau1 [s], else from the mean induction and the wind speed
    };

    double oye_tau1(double a_mean, double R, double U0);
    double oye_tau2(double tau1, double r_over_R);

    struct inflow_state
    {
        double Wqs[2] = {0.0, 0.0};     // quasi-steady induced velocity at the committed time (normal, tangential)
        double Wint[2] = {0.0, 0.0};
        double W[2] = {0.0, 0.0};       // the dynamic induced velocity
    };

    // steady: all equal to W_qs
    inflow_state inflow_steady(const double Wqs[2]);
    // the state at t0 + dt from the state at t0 and W_qs at t0 + dt
    inflow_state inflow_step(const inflow_state&, const double Wqs[2], double tau1, double tau2, double k, double dt);

    struct dynamic_stall_settings
    {
        int model = 0;                  // 0 off, 1 Øye, 2 Beddoes-Leishman (HGM)
        // < 0: the airfoil file's value, else the default above (T_f0: 6 with Øye, 3 with HGM)
        double A1 = -1.0, A2 = -1.0, b1 = -1.0, b2 = -1.0, T_f0 = -1.0, T_p = -1.0;
    };

    // the static polar prepared for dynamic stall (one airfoil, the table at the Reynolds number re)
    class stall_polar
    {
    public:
        stall_polar() = default;
        stall_polar(const airfoil&, double re, const dynamic_stall_settings&);

        double alpha0 = 0.0, cl_alpha = 2.0*3.14159265358979323846;
        double A1 = 0.165, A2 = 0.335, b1 = 0.0455, b2 = 0.3, T_f0 = 3.0, T_p = 1.7;
        airfoil_coef stat(double alpha) const;          // the static coefficients
        double f_st(double alpha) const;
        double cl_fs(double alpha) const;
        double cl_inv(double alpha) const {return cl_alpha*(alpha - alpha0);}
        double cd0() const {return cd0_;}

    private:
        const airfoil *af_ = nullptr;
        double re_ = 0.0, cd0_ = 0.0;
        std::vector<double> a_, f_, fs_;                // f_st and Cl_fs at the table's angles (rad, increasing)
        double interp(const std::vector<double>&, double alpha) const;
    };

    struct stall_input
    {
        double alpha = 0.0;     // at the aerodynamic centre [rad]
        double U = 0.0;         // relative speed [m/s]
        double chord = 0.0;
        double omega = 0.0;     // nose-up angular velocity of the section [rad/s]
    };

    struct stall_state
    {
        double x[4] = {0.0, 0.0, 0.0, 1.0};   // Øye: x[3] = f
        double a34 = 0.0, clp = 0.0, fF = 1.0, fst = 1.0;   // the inputs at the committed time
    };

    double stall_Tu(const stall_input&);
    // steady state at the input (the coefficients then are the static ones)
    stall_state stall_steady(const stall_polar&, int model, const stall_input&);
    // the state at t0 + dt from the state at t0 and the input at t0 + dt; the coefficients at t0 + dt
    airfoil_coef stall_step(const stall_polar&, int model, const stall_state&, const stall_input&, double dt, stall_state &out);
}

#endif
