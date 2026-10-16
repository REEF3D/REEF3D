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
// axis (revolute, cylindrical) between the two bodies of the joint. With q the joint coordinate
// (zero at creation) and qd its rate:
//   F = F0 - k (q - q0) - c qd - Fc tanh(qd / v_eps)      (+ end stops)
// limited to |F| <= F_max and |F qd| <= P_max; absorbed power P = -F qd.
// End stops: penalty spring k_stop and damper c_stop beyond q_min / q_max.
// Control: passive; latching (held with c_latch for t_hold after each zero crossing of qd);
// declutching (free, F = 0 apart from the stops, for t_hold after each zero crossing of qd).
// Positive F pushes body b along the axis relative to body a (increases q).

#ifndef CORAL_PTO_H_
#define CORAL_PTO_H_

#include"coral/core/coral_load.h"
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
        double c_latch = 0.0;
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
        double power(const model &m) const {return -force(m)*rate(m);}    // absorbed power
        double energy() const {return energy_;}            // absorbed energy so far
        double mean_power(double t) const {return t > t0_ ? energy_/(t - t0_) : 0.0;}
        bool held() const {return held_;}

        pto_settings settings;

    private:
        std::string name_;
        const joint *j_;
        bool rot_;
        bool held_ = false;          // latched or declutched
        bool armed_ = true;          // a zero crossing may start the next hold
        double t_release_ = 0.0;
        double qd_last_ = 0.0;
        double P_last_ = 0.0;
        double energy_ = 0.0;
        double t0_ = 0.0;
    };
}

#endif
