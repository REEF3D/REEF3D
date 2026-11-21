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
// REEF3D::CORAL: hydraulic power take-off (PTO-Sim type): a double-acting cylinder in the PTO joint
// (or link) rectified by four check valves into a high-pressure (HP) and out of a low-pressure (LP)
// gas accumulator, a hydraulic motor between them that drives the generator shaft, and a relief valve
// from HP to LP. Pressures relative to the ambient, x = q the piston stroke (zero at the start),
// v = qd. Positive F pushes body b along +q (extends the cylinder).
//
//   chambers      V_a = V_a0 + A_a x,  V_b = V_b0 - A_b x
//                 p_a' = beta/V_a (-A_a v - Q_a,hp + Q_lp,a),  p_b' = beta/V_b (A_b v - Q_b,hp + Q_lp,b)
//                 F = A_a p_a - A_b p_b
//   check valves  Q = Cd A_v s((dp - dp_crack)/dp_open) sqrt(2 |dp|/rho) for dp > dp_crack,
//                 s(u) = 3u^2 - 2u^3 (0..1, smooth opening), dp across the valve in its flow direction
//   accumulators  p = p_pre (V/V_g)^kappa (gas volume V_g, total volume V, precharge p_pre),
//                 V_g,hp' = -(Q_a,hp + Q_b,hp - Q_m - Q_r),  V_g,lp' = -(Q_m + Q_r - Q_lp,a - Q_lp,b)
//   motor         Q_m = D w + C_leak (p_hp - p_lp),  T_m = eta_m D (p_hp - p_lp) (1/eta_m when it pumps)
//   relief        Q_r = K_r max(0, p_hp - p_relief)
//   shaft         J w' = T_m - T_gen - c_w w,  T_gen = b_gen w (|T_gen| <= T_max if T_max > 0)
// The generator torque law b_gen is what a controller sets (set_param "b_gen"); the generator model
// (coral_generator.h) turns T_gen w into electrical power.
//
// The valves make the system stiff (the chambers settle on the valve flow in microseconds), so its
// states are integrated implicitly (two-stage SDIRK, L-stable and second order, Newton in each stage,
// sub-steps dt_sub) from the state of
// the last accepted step to the time of each force() call along the cubic Hermite path of the piston
// between the accepted (q, qd) and the asked one. The force is thus a function of the state of the
// step (implicit coupling with the multibody Newton iteration); accept() moves the reference state on.

#ifndef CORAL_PTO_HYDRAULIC_H_
#define CORAL_PTO_HYDRAULIC_H_

#include"coral/pto/coral_pto_element.h"
#include"coral/pto/coral_generator.h"
#include<array>

namespace coral
{
    struct hydraulic_settings
    {
        double A_a = 0.0, A_b = 0.0;            // piston areas, cap and rod side [m^2]
        double V_a0 = 0.0, V_b0 = 0.0;          // chamber volumes at x = 0 with the hoses [m^3]
        double beta = 1.5e9;                    // effective bulk modulus [Pa]
        double rho = 850.0;                     // oil density [kg/m^3]
        double Cd = 0.61, A_v = 0.0;            // check valves: discharge coefficient, open area [m^2]
        double dp_crack = 0.0, dp_open = 1.0e5; // cracking pressure, opening range [Pa]
        double V_hp = 0.0, p_hp_pre = 0.0, p_hp0 = 0.0;   // HP accumulator: volume, precharge, initial pressure
        double V_lp = 0.0, p_lp_pre = 0.0, p_lp0 = 0.0;   // LP accumulator
        double kappa = 1.4;                     // polytropic exponent of the gas
        double D = 0.0;                         // motor displacement [m^3/rad]
        double C_leak = 0.0;                    // motor leakage [m^3/(s Pa)]
        double eta_m = 1.0;                     // motor torque efficiency
        double J = 1.0, c_w = 0.0;              // shaft inertia [kg m^2], friction [N m s]
        double w0 = 0.0;                        // initial shaft speed [rad/s]
        double b_gen = 0.0, T_max = 0.0;        // generator torque law T = b_gen w, its limit (0: none)
        double p_relief = 0.0, K_r = 0.0;       // relief valve HP -> LP (K_r = 0: none) [m^3/(s Pa)]
        double dt_sub = 1.0e-3;                 // sub-step of the implicit integration [s]
        bool has_gen = false;                   // electrical generator on the shaft
        generator_settings gen;

        // empty if fine, else the problem
        std::string check() const;
    };

    class pto_hydraulic : public pto_base
    {
    public:
        // states: p_a, p_b, V_g,hp, V_g,lp, w
        using state = std::array<double,5>;

        explicit pto_hydraulic(const hydraulic_settings&);

        pto_output force(const pto_state&) override;
        void accept(const pto_state&, double) override;
        void start(const pto_state&) override;
        void set_param(const std::string&, double) override;
        const char* name() const override {return "hydraulic";}
        bool scaled() const override {return false;}
        void channel_names(const std::string&, std::vector<std::string>&) const override;
        void channel_values(std::vector<double>&) const override;
        void save_state(std::vector<double>&) const override;
        void load_state(const std::vector<double>&, size_t&) override;

        // pressures of the accumulators for gas volumes
        double p_hp(double Vg) const;
        double p_lp(double Vg) const;
        // the state of the last force() call and the accepted one
        const state &last() const {return y_;}
        const state &accepted() const {return yc_;}
        // energies of the accepted steps: piston work in, generator shaft work out, dissipated in valves,
        // leakage, motor, shaft friction and relief, electrical energy delivered
        double E_in = 0.0, E_gen = 0.0, E_diss = 0.0, E_el = 0.0;
        // energy stored in the gas of both accumulators relative to the start (polytropic work)
        double gas_energy(const state&) const;
        double shaft_energy(const state &y) const {return 0.5*set.J*y[4]*y[4];}

        hydraulic_settings set;

    private:
        struct flows {double Qa_hp, Qb_hp, Qlp_a, Qlp_b, Qm, Qr, Tm, Tg, P_diss;};
        flows eval_flows(double x, const state&) const;
        void rhs(double x, double v, const state&, state &dy, flows *fl=nullptr) const;
        double valve(double dp) const;
        double gen_torque(double w) const;
        // integrates y from the accepted state to s; energies of the path when e != nullptr
        void integrate(const pto_state &s, state &y, double *e) const;
        double piston_force(const state &y) const {return set.A_a*y[0] - set.A_b*y[1];}

        bool started_ = false;
        state yc_{}, y_{};
        double tc_ = 0.0, qc_ = 0.0, qdc_ = 0.0;     // accepted time, stroke and rate
        double Pel_last_ = 0.0;
        // cache of the last force() call
        double t_last_ = -1.0e300, q_last_ = 0.0, qd_last_ = 0.0;
        pto_output out_last_;
        generator_model gen_;
    };
}

#endif
