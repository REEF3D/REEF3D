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
// REEF3D::CORAL: power take-off in a joint: axial force (prismatic, cylindrical) or torque about the
// axis (revolute, cylindrical) between the two bodies of the joint, or between two points (link: a
// point on body a and one on body b, q = L - L0 the change of their distance, the force along the
// line between them, positive pushing them apart; an actuator or hydraulic cylinder mounted between
// a flap and the ground, a rod, a tether to a winch: the motion conversion of the geometry is
// exact), built from the elements of the PTO library (coral_pto_element.h). With q the joint
// coordinate (zero at creation) and qd its rate:
//   generator  F = F0 - k (q - q0) - c qd - Fc tanh(qd / v_eps),  |F| <= F_max, |F qd| <= P_max
//   losses     F = - c_loss qd - Fc_loss tanh(qd / v_eps)   (drive train: belt, bearings, seals)
//   end stops  penalty spring k_stop and damper c_stop beyond q_min / q_max (push only)
// absorbed power P = -F qd of the generator (the end stops and the latch absorb nothing); the
// losses act on the bodies in every state (also latched or declutched, not limited) and their power
// is dissipated (P_loss, E_loss), not absorbed.
// Control: passive; latching (a latch element, spring k_latch and damper c_latch at the position of
// the turn, holds the joint for t_hold after each zero crossing of qd); declutching (generator off,
// F = 0 apart from the stops, for t_hold after each zero crossing of qd). The crossing times are
// interpolated between the steps. The controllers act once per converged step (commit).
// Positive F pushes body b along the axis relative to body a (increases q).
// Hydraulic PTO (coral_pto_hydraulic.h): cylinder, check valves, accumulators, motor and generator
// shaft as an element of its own (not scaled or limited with the generator, its power absorbed).
// Generator (coral_generator.h): the electrical power of the generator force (direct drive, through
// a gear ratio) or of the shaft of the hydraulic motor; channels _Pe, _Ee.
// Further controllers (coral_pto_wec_control.h, add_controller): gain schedule, user library, model
// predictive control; they set the gains k, c and the force F0 of the generator once per step.

#ifndef CORAL_PTO_H_
#define CORAL_PTO_H_

#include"coral/core/coral_load.h"
#include"coral/pto/coral_pto_element.h"
#include"coral/pto/coral_pto_controller.h"
#include"coral/pto/coral_pto_hydraulic.h"
#include"coral/pto/coral_generator.h"
#include<memory>
#include<limits>
#include<string>
#include<vector>

namespace coral
{
    class joint;

    enum class pto_control {passive, latching, declutching};

    struct pto_settings
    {
        double k = 0.0, c = 0.0, q0 = 0.0, F0 = 0.0;
        double coulomb = 0.0, v_eps = 1.0e-3;
        double c_loss = 0.0, coulomb_loss = 0.0;    // drive-train losses, 0: none
        double F_max = 0.0;          // 0: no force limit
        double P_max = 0.0;          // 0: no power limit
        double q_min = -std::numeric_limits<double>::infinity();
        double q_max = std::numeric_limits<double>::infinity();
        double k_stop = 0.0, c_stop = 0.0;
        pto_control control = pto_control::passive;
        double t_hold = 0.0;
        double c_latch = 0.0, k_latch = 0.0;
        bool commanded = false;      // a controller sets k, c, F0: the spring-damper is there even if all are 0
        bool has_hydraulic = false;  // hydraulic PTO element (its generator: hydraulic.gen)
        hydraulic_settings hydraulic;
        bool has_generator = false;  // generator of the direct drive (without the hydraulic PTO)
        generator_settings generator;
    };

    class pto_load : public load
    {
    public:
        // the joint must live as long as the load (it is owned by the model)
        pto_load(const joint&, bool rotational, const pto_settings&, const std::string &name="pto");
        // link between the global points pa on body a and pb on body b (-1: ground), in the configuration
        // when the PTO is created; q = |pb - pa| - L0, L0 the distance at creation
        pto_load(const model&, int a, const Vec3 &pa, int b, const Vec3 &pb, const pto_settings&,
                 const std::string &name="pto");

        // a further controller, acting once per converged step after the built-in ones (kept by initialize())
        void add_controller(std::unique_ptr<pto_controller> c) {ext_.push_back(std::move(c));}

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        // restart from a checkpoint (coral_state.h)
        void save_state(std::vector<double>&) const override;
        bool load_state(const std::vector<double>&, size_t &k) override;
        double potential(const model&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        double coordinate(const model&) const;
        double rate(const model&) const;
        double force(const model&) const;      // force (torque) of the PTO in the current state
        double force(const model&, double t) const;
        double power(const model &m) const;    // absorbed power of the generator
        double loss_power(const model&) const;   // power dissipated in the drive-train losses
        double energy() const {return energy_;}            // absorbed energy so far
        double loss_energy() const {return energy_loss_;}  // dissipated energy so far
        bool has_loss() const {return settings.c_loss > 0.0 || settings.coulomb_loss > 0.0;}
        bool has_generator() const {return settings.has_generator && !settings.has_hydraulic;}
        double electrical_energy() const;      // delivered electrical energy (direct or hydraulic generator)
        // electrical power of the direct-drive generator in the current state (0 without one)
        generator_power electrical(const model&) const;
        // the hydraulic element (nullptr without one)
        const pto_hydraulic *hydraulic() const {return hyd_;}
        bool is_link() const {return j_ == nullptr;}
        // apply the force F along the joint or link to the generalized forces
        void apply(const model&, double F, VecX &f) const;
        double mean_power(double t) const {return t > t0_ ? energy_/(t - t0_) : 0.0;}
        bool held() const;                     // latched or declutched

        // the elements (built from the settings in initialize())
        const pto_composite &composite() const {return pto_;}

        pto_settings settings;

    private:
        void build();
        pto_output evaluate(const model&, double t) const;

        void update_controllers(const model&, double t, double q, double qd);

        std::string name_;
        const joint *j_ = nullptr;
        bool rot_ = false;
        int a_ = -1, b_ = -1;        // link: bodies, points in their frames, initial distance
        Vec3 sa_ = Vec3::Zero(), sb_ = Vec3::Zero();
        double L0_ = 0.0;
        mutable pto_composite pto_;
        std::vector<std::unique_ptr<pto_controller>> ctrl_;
        std::vector<std::unique_ptr<pto_controller>> ext_;
        pto_hydraulic *hyd_ = nullptr;
        generator_model gen_;
        double Pe_last_ = 0.0, E_el_ = 0.0;
        double t_ = 0.0;             // time of the last converged state
        double P_last_ = 0.0;
        double energy_ = 0.0;
        double Ploss_last_ = 0.0;
        double energy_loss_ = 0.0;
        double t0_ = 0.0;
    };
}

#endif
