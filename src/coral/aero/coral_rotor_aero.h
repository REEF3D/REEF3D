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
// momentum at each aerodynamic station from the station's actual motion.
//
// A station sits on a body (the rotor, or later a blade or a beam node) at a point with three directions fixed to
// it: n (normal to the span in the plane through the shaft axis, downwind), t (the direction the blade moves in when
// the rotor turns forwards) and the span direction. At time t:
//   V  = u_wind(P, t) + u_tower(P) - v_P                       (relative inflow; the tower's potential flow optional)
//   Vx = V.n,  Vy = -V.t,  r = distance of P from the shaft axis,  theta = twist + pitch
//   bem_solve (coral_bem.h) with R and R_hub in the rotor plane (cone included)
//   skewed wake (Pitt-Peters, as AeroDyn's SkewMod 1): gamma = angle between the shaft axis and the inflow at the hub,
//   chi = (0.6 a + 1) gamma, the induction a (1 + 15 pi/32 (r/R) tan(chi/2) cos psi), psi the azimuth of the station
//   from the downwind side of the cross flow
//   loads: F = w (fn n + ft t), M = w mz (n x t) at P, w the station's share of the span (trapezoid rule)
// Tower potential flow (AeroDyn's TwrPotent 1): around a tower from base to top with linear diameter D, the wind
// normal to the tower axis U, x along it and y across, scaled by D/2: u = U (1 - (x^2 - y^2)/(x^2 + y^2)^2),
// v = -2 U x y/(x^2 + y^2)^2 outside the tower, between its base and top.
// Rotor values: thrust T = sum F.e, torque Q = sum ((P - H) x F + M).e about the hub centre H, power Q omega, Cp and Ct
// with the wind speed at the hub and the swept area in the rotor plane; aerodynamic root moments of each blade about
// the root, out of the rotor plane (about span x n, positive for downwind thrust) and in it (about span x t, positive
// for a driving force).

#ifndef CORAL_ROTOR_AERO_H_
#define CORAL_ROTOR_AERO_H_

#include"coral/core/coral_load.h"
#include"coral/aero/coral_bem.h"
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
        Vec3 n0 = Vec3::UnitX(), t0 = Vec3::UnitY(), z0 = Vec3::UnitZ();   // directions in the body frame
        double chord = 0.0, twist = 0.0, w = 0.0;
        int airfoil = 0, blade = 0;
    };

    struct rotor_aero_settings
    {
        int B = 3;
        double R_hub = 0.0, R_tip = 1.0;    // in the rotor plane
        double rho = 1.225;
        bem_options bem;
        bool skew = true;
        // tower potential flow (on the tower body, -1: fixed)
        bool tower = false;
        int tower_body = -1;
        Vec3 tower_base = Vec3::Zero(), tower_top = Vec3::UnitZ();     // global at the set-up
        double D_base = 0.0, D_top = 0.0;
    };

    struct station_load
    {
        Vec3 P, F, M;
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
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

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

    private:
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
    };
}

#endif
