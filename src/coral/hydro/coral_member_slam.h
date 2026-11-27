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
// REEF3D::CORAL: slamming of breaking waves on a flexible member (a monopile), in addition to its Morison or
// slender-body loads (coral_member_morison.h). The line force of an impact with velocity V on the cylinder of radius R
//   f(t) = rho R V^2 Cs(tau),   tau = V (t - t0) / R,
// acts uniformly over the curling part of the crest, from (1 - lambda) eta_b to eta_b above the still water level
// (eta_b: the crest elevation at the member), in the horizontal direction of the crest velocity, normal to the member.
// Slamming coefficients (model):
//   wo      Wienke and Oumeraci (2005), Coastal Engineering 52, 435-462, as printed in Ma et al. (2023), J. Marine Sci.
//           Appl. 22, 1-13, Eq. 3:
//             Cs = 2 pi - 2 sqrt(tau) artanh sqrt(1 - tau/4),                                    0 <= tau <= 1/8
//             Cs = pi / sqrt(6 tau') - (8 tau'/3)^(1/4) artanh sqrt(1 - tau' sqrt(6 tau')),    tau' = tau - 1/32,
//                                                                                              1/8 < tau <= 13/32
//           (the published form jumps at tau = 1/8, 4.57 -> 2.77, and ends at 1.30 at tau = 13/32)
//   goda    Goda et al. (1966), von Karman: Cs = pi (1 - tau),  0 <= tau <= 1
//   cw      Campbell and Weynberg (1980), DNV-RP-C205: Cs = 5.15 (D/(D + 19 s) + 0.107 s/D), penetration s = V (t - t0),
//           until the member is immersed, s = D (tau = 2)
//   custom  Cs = cs0 (1 - tau/tau_end) (shape linear) or cs0 (shape constant) for 0 <= tau <= tau_end: a fitted model
//           (e.g. a slamming coefficient and duration from model tests)
// Events (trigger), evaluated once per converged step (commit) at the member's waterline point:
//   breaking   the resolved fluid flags breaking at the member (REEF3D::FNPF's breaking model, A 350) while the crest
//              rises; the impact starts once the crest has passed the member (eta falls), at the crest's elevation
//              eta_b and velocity
//   kinematic  the same with the kinematic criterion B = |u_h| / c >= Bc at the surface (Barthelemy et al. 2018:
//              Bc = 0.85), c from the linear dispersion relation of period T at the depth
//   time       impacts at given times (eta_b and V given, or taken at that time; without waves given, along x)
// One impact per crest: the trigger is armed again once the surface falls below the still water level.
// Impact velocity V: the celerity c of the period T at the depth (default; Wienke and Oumeraci, IEC 61400-3: the breaking
// wave's celerity), vfac |u_h| (velocity crest: u_h the horizontal fluid velocity at the crest, the surface velocity of a
// sampled fluid) or given; the direction is that of u_h.
// The impact starts one step after the crest passed the member's axis, i.e. about R/c + dt after the first contact.

#ifndef CORAL_MEMBER_SLAM_H_
#define CORAL_MEMBER_SLAM_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_points.h"
#include"coral/hydro/coral_waves.h"
#include<memory>
#include<string>
#include<utility>
#include<vector>

namespace coral
{
    enum class slam_model {wo, goda, cw, custom};
    enum class slam_trigger {breaking, kinematic, time};

    struct member_slam_coef
    {
        double D = 1.0;
        slam_model model = slam_model::wo;
        double lambda = 0.46;          // curling factor
        double cs0 = 2.0*3.14159265358979323846;    // custom: Cs at the start
        double tau_end = 1.0;          // custom: duration in R/V
        bool constant = false;         // custom: constant Cs (default linear decay to 0)
        slam_trigger trigger = slam_trigger::breaking;
        double Bc = 0.85;              // kinematic: breaking index
        double period = -1.0;          // period of the celerity c (kinematic trigger, impact velocity)
        double depth = -1.0;           // depth of the celerity (< 0: the sampled bed or the waves' depth)
        std::vector<double> times;     // time: impact times
        double V = -1.0;               // impact velocity (< 0: vfac c, or vfac |u_h| at the crest)
        bool crest_velocity = false;
        double vfac = 1.0;
        double eta_b = -1.0;           // time: crest elevation (< 0: the elevation at that time)
        int n = 4;                     // sub-strips of the impact zone per element
    };

    // the slamming coefficient and the duration (in tau = V t / R) of a model
    double slam_coefficient(const member_slam_coef&, double tau);
    double slam_duration(const member_slam_coef&);

    class member_slam : public load
    {
    public:
        member_slam(std::shared_ptr<const point_set> pts, const std::vector<std::pair<int,int>> &elements,
                    const member_slam_coef&, std::shared_ptr<const wave_kinematics>, double rho, double g, const std::string &name);

        std::string name() const override {return name_;}
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.insert(b.end(), bodies_.begin(), bodies_.end());}
        void coupled_coordinates(std::vector<int> &c) const override {c.insert(c.end(), coords_.begin(), coords_.end());}
        // the impact load depends on time only (its weak dependence on the member's position is left out of the Jacobian)
        bool has_tangent() const override {return true;}
        void kinematics_points(const model&, const void *source, std::vector<Vec3>&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;

        struct event {double t0, V, eta_b; Vec3 dir;};
        const std::vector<event> &events() const {return ev_;}
        int nevents() const {return nev_;}
        // total force and moment about the origin at time t
        void total(const model&, double t, Vec3 &F, Vec3 &M) const;

    private:
        template<class Sink> void impacts(const model&, double t, Sink sink) const;
        Vec3 waterline(const model&) const;
        void find_waterline(const model&) const;
        void fire(double t0, double eta_b, const Vec3 &uh, double cel);
        double celerity(const Vec3 &p) const;
        std::string name_;
        std::shared_ptr<const point_set> pts_;
        std::vector<std::pair<int,int>> elem_;
        std::vector<int> pidx_;
        member_slam_coef c_;
        std::shared_ptr<const wave_kinematics> waves_;
        double rho_, g_;
        std::vector<int> bodies_, coords_;
        mutable int iwl_ = -1;         // element crossing z = 0 at the start
        // trigger state
        std::vector<event> ev_;
        int nev_ = 0;
        bool armed_ = true, seen_ = false;
        int nt_ = 0;                   // next impact time (time)
        double eta1_ = 0.0, eta2_ = 0.0, B_ = 0.0, eta_ = 0.0;
        Vec3 uh1_ = Vec3::Zero();
        int nstep_ = 0;
    };
}

#endif
