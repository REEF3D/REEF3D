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
// REEF3D::CORAL: aerodynamic load on the blades of a rotor in the multibody model, quasi-steady blade element
// momentum at each aerodynamic station from the station's actual motion and orientation.
//
// A station sits on a body (the rotor, a pitch bearing or a beam node of a flexible blade) at a point, with the span
// direction z and the chord direction c (leading edge, at zero pitch) fixed to that body. At time t, with the shaft
// axis e (downwind) and the sense of rotation s (+1: positive about e):
//   n = (e - (e.z) z)/|.|  normal to the span in the plane through the axis (downwind),  t = s (e x z)/|.|  the
//   direction of motion of the blade;  theta = atan2(-c.n, c.t) + pitch  (the chord's angle from the rotor plane,
//   twist, pitch of a pitch joint and elastic torsion included; pitch: the collective pitch of rigid blades without
//   pitch joints)
//   V  = u_wind(P, t) + u_tower(P) - v_P,  Vx = V.n,  Vy = -V.t,  r = distance of P from the shaft axis
//   bem_solve (coral_bem.h) with R and R_hub in the rotor plane (cone included)
//   skewed wake (Pitt-Peters, as AeroDyn's SkewMod 1): gamma = angle between the shaft axis and the inflow at the hub,
//   chi = (0.6 a + 1) gamma, the induction a (1 + 15 pi/32 (r/R) tan(chi/2) cos psi), psi the azimuth of the station
//   from the downwind side of the cross flow
//   loads: F = w (fn n + ft t), M = w mz (t x n) at P (mz nose-up positive), w the station's share of the span
// A station at the tip or at the hub of the set-up geometry keeps r = R_tip or R_hub (no load with the tip or hub
// loss), whatever rounding in the turned frames or the blade's deflection does to its radius.
// Tower potential flow (AeroDyn's TwrPotent 1): around a tower from base to top with linear diameter D, the wind
// normal to the tower axis U, x along it and y across, scaled by D/2: u = U (1 - (x^2 - y^2)/(x^2 + y^2)^2),
// v = -2 U x y/(x^2 + y^2)^2 outside the tower, between its base and top.
// Tangent: the aerodynamic damping, each station's generalized force against the velocity and angular velocity of its
// own body, from the derivative of its force by its relative inflow (forward differences: 4 station evaluations), and
// against the rotation of its body about the span (pitch, torsion: 2 evaluations), not of the whole rotor. The dependence on positions and on
// the hub and the shaft axis is left out (an approximate Jacobian: Newton converges to the same residual). With the
// position part, a blade turned out of the plane of the (fixed) shaft axis gave a spurious stiffness and Newton needed
// twice the iterations on a flexible tower; without the whole rotor in the finite-difference Jacobian a turbine on a
// flexible monopile runs 3.5 times faster.
// Unsteady aerodynamics (coral_unsteady_aero.h; settings dyn_inflow, dyn_stall): after the quasi-steady solution of a
// station, its induced velocities W_qs = (a Vx, a' Vy) (skewed wake included) pass Øye's dynamic inflow filter, and the
// angle of attack from the filtered inflow (Vx - W_n, Vy + W_t) passes the dynamic stall model (Øye or Beddoes-Leishman
// HGM, with the section's pitch and torsion rate relative to the hub and the rate of a pitch table); the loads follow
// from the resulting cl, cd, cm (as AeroDyn: the momentum balance with the static polar, the unsteady models on the
// loads). The states are committed once per converged step (commit) and each trial state of a step is integrated from
// them (eval, tangent: the aerodynamic damping then sees the lag of the induction and of the separation); tau1 from
// the area-weighted mean axial induction and the wind at the hub at the committed state. Before initialize() (static
// equilibrium, modes) and without the models: quasi-steady. save_state/load_state carry the states for a restart.
// Rotor values: thrust T = sum F.e, torque Q = sum ((P - H) x F + M).e about the hub centre H, power Q omega, Cp and Ct
// with the wind speed at the hub and the swept area in the rotor plane; aerodynamic root moments of each blade about
// the root, out of the rotor plane (about span x n, positive for downwind thrust) and in it (about span x t, positive
// for a driving force).

#ifndef CORAL_ROTOR_AERO_H_
#define CORAL_ROTOR_AERO_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_body.h"
#include"coral/aero/coral_bem.h"
#include"coral/aero/coral_unsteady_aero.h"
#include"coral/wind/coral_wind.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    class joint;

    struct aero_station
    {
        int body = -1;
        Vec3 s = Vec3::Zero();                                  // point in the body frame
        Vec3 c0 = Vec3::UnitY(), z0 = Vec3::UnitZ();            // chord (leading edge) and span in the body frame
        double chord = 0.0, w = 0.0;
        int airfoil = 0, blade = 0;
    };

    struct rotor_aero_settings
    {
        int B = 3;
        double sense = 1.0;                 // +1: the rotor turns positively about the shaft axis
        double R_hub = 0.0, R_tip = 1.0;    // in the rotor plane
        double rho = 1.225;
        bem_options bem;
        bool skew = true;
        // tower potential flow (on the tower body, -1: fixed)
        bool tower = false;
        int tower_body = -1;
        Vec3 tower_base = Vec3::Zero(), tower_top = Vec3::UnitZ();     // global at the set-up
        double D_base = 0.0, D_top = 0.0;
        // unsteady aerodynamics (coral_unsteady_aero.h)
        dynamic_inflow_settings dyn_inflow;
        dynamic_stall_settings dyn_stall;
    };

    struct station_load
    {
        Vec3 P, F, M;
        Vec3 n, t;          // rotor-plane directions at the station
        double theta = 0.0;
        bem_result bem;
    };

    class rotor_aero : public load
    {
    public:
        // hub: the rotor body; H the hub centre (global, at the set-up); shaft: its joint (axis along the downwind
        // shaft); roots: the root point of each blade (global)
        rotor_aero(const model&, std::shared_ptr<const wind_inflow>, std::vector<airfoil>, int hub, const Vec3 &H,
                   const joint &shaft, std::vector<aero_station>, const std::vector<Vec3> &roots,
                   const rotor_aero_settings&, const std::string &name="rotor");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        bool has_tangent() const override {return true;}
        void tangent(const model&, double t, double cK, double cC, double cM, std::vector<triplet>&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;

        // pitch of each blade [rad], or a collective table in time (rad)
        void set_pitch(int blade, double p) {pitch_[blade] = p;}
        void set_pitch_table(std::vector<double> t, std::vector<double> p) {pt_ = std::move(t); pv_ = std::move(p);}
        double pitch(int blade, double t) const;

        // the loads of all stations at the current state
        void station_loads(const model&, double t, std::vector<station_load>&) const;
        struct rotor_values {double T, Q, P, Cp, Ct, lambda, U, yaw, omega, azimuth;};
        rotor_values values(const model&, double t, const std::vector<station_load>&) const;
        // the velocity the tower's potential flow adds to the wind u at P
        Vec3 tower_influence(const model&, const Vec3 &P, const Vec3 &u) const;

        Vec3 hub_centre(const model&) const;
        const std::vector<aero_station> &stations() const {return st_;}
        const rotor_aero_settings &settings() const {return set_;}
        bool unsteady() const {return set_.dyn_inflow.on || set_.dyn_stall.model > 0;}
        double tau1() const {return tau1_;}
        // the committed unsteady states of the stations (in the order of stations())
        const std::vector<inflow_state> &inflow_states() const {return di_;}
        const std::vector<stall_state> &stall_states() const {return ds_;}

    private:
        // the inflow state at the hub shared by all stations
        struct hub_state {Vec3 H, e; double gam; Vec3 cross;};
        hub_state hub(const model&, double t) const;
        // one station on the body b (its state as given); i >= 0: the station's index, for the unsteady aerodynamics
        // (with di, ds: their states at t)
        void station(const model&, const rigid_body &b, const aero_station&, const hub_state&, double t, station_load&,
                     const Vec3 &dV=Vec3::Zero(), int i=-1, inflow_state *di=nullptr, stall_state *ds=nullptr) const;
        // the unsteady aerodynamics of a station after its quasi-steady solution (el, Vx, Vy: its element and inflow;
        // c: its chord direction, global)
        void unsteady_station(const model&, const rigid_body &b, const aero_station&, const hub_state&, double t,
                              station_load&, const bem_element &el, double Vx, double Vy, const Vec3 &c, int i,
                              inflow_state *di, stall_state *ds) const;
        // the unsteady states at t of all stations: steady (initialize) or a step from the committed ones
        void advance(const model&, double t);

        std::string name_;
        std::shared_ptr<const wind_inflow> wind_;
        std::vector<airfoil> af_;
        int hub_;
        Vec3 sH_;
        const joint *shaft_;
        std::vector<aero_station> st_;
        std::vector<Vec3> roots_;          // in the hub frame
        rotor_aero_settings set_;
        Vec3 tb_, tt_;                     // tower in its body frame
        std::vector<double> pitch_, pt_, pv_;
        std::vector<char> edge_;           // per station: 1 at the tip, 2 at the hub (set-up geometry), else 0
        // unsteady aerodynamics: the committed states at tc_ (started_ after initialize; before it, and without the
        // unsteady models, quasi-steady)
        bool started_ = false;
        double tc_ = 0.0, tau1_ = 0.0;
        std::vector<inflow_state> di_;
        std::vector<stall_state> ds_;
        std::vector<stall_polar> sp_;     // per airfoil
    };
}

#endif
