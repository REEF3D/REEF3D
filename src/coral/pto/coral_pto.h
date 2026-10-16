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
// axis (revolute, cylindrical) between the two bodies of the joint, built from the elements of the
// PTO library (coral_pto_element.h). With q the joint coordinate (zero at creation) and qd its rate:
//   generator  F = F0 - k (q - q0) - c qd - Fc tanh(qd / v_eps),  |F| <= F_max, |F qd| <= P_max
//   end stops  penalty spring k_stop and damper c_stop beyond q_min / q_max (push only)
// absorbed power P = -F qd of the generator (the end stops and the latch absorb nothing).
// Control: passive; latching (a latch element, spring k_latch and damper c_latch at the position of
// the turn, holds the joint for t_hold after each zero crossing of qd); declutching (generator off,
// F = 0 apart from the stops, for t_hold after each zero crossing of qd). The crossing times are
// interpolated between the steps. The controllers act once per converged step (commit).
// Positive F pushes body b along the axis relative to body a (increases q).

#ifndef CORAL_PTO_H_
#define CORAL_PTO_H_

#include"coral/core/coral_load.h"
#include"coral/pto/coral_pto_element.h"
#include"coral/pto/coral_pto_controller.h"
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
        double F_max = 0.0;          // 0: no force limit
        double P_max = 0.0;          // 0: no power limit
        double q_min = -std::numeric_limits<double>::infinity();
        double q_max = std::numeric_limits<double>::infinity();
        double k_stop = 0.0, c_stop = 0.0;
        pto_control control = pto_control::passive;
        double t_hold = 0.0;
        double c_latch = 0.0, k_latch = 0.0;
    };

    class pto_load : public load
    {
    public:
        // the joint must live as long as the load (it is owned by the model)
        pto_load(const joint&, bool rotational, const pto_settings&, const std::string &name="pto");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void initialize(const model&, double t0) override;
        void eval(const model&, double t, VecX &f) const override;
        void commit(const model&, double t, double dt) override;
        double potential(const model&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        double coordinate(const model&) const;
        double rate(const model&) const;
        double force(const model&) const;      // force (torque) of the PTO in the current state
        double force(const model&, double t) const;
        double power(const model &m) const;    // absorbed power of the generator
        double energy() const {return energy_;}            // absorbed energy so far
        double mean_power(double t) const {return t > t0_ ? energy_/(t - t0_) : 0.0;}
        bool held() const;                     // latched or declutched

        // the elements (built from the settings in initialize())
        const pto_composite &composite() const {return pto_;}

        pto_settings settings;

    private:
        void build();
        pto_output evaluate(const model&, double t) const;

        std::string name_;
        const joint *j_;
        bool rot_;
        mutable pto_composite pto_;
        std::vector<std::unique_ptr<pto_controller>> ctrl_;
        double t_ = 0.0;             // time of the last converged state
        double P_last_ = 0.0;
        double energy_ = 0.0;
        double t0_ = 0.0;
    };
}

#endif
