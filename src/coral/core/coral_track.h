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
// REEF3D::CORAL: the motion of bodies over a time step, recorded at the ends of the (sub)steps, and
// interpolated in between for a coupled solver that needs the bodies at its own times (the stages of
// the fluid's Runge-Kutta scheme). Cubic Hermite interpolation per interval [t_a, t_b]:
//   position from x and v, velocity from v and the acceleration, rotation R(t) = exp(theta(t)) R_a with
//   theta from 0, log(R_b R_a^T) and the global angular velocities at both ends, angular velocity from
//   omega and the angular acceleration (all global). At the recorded times the states are exact.

#ifndef CORAL_TRACK_H_
#define CORAL_TRACK_H_

#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    class model;

    struct body_snapshot
    {
        Vec3 x, v, a;           // CoG: position, velocity, acceleration (global)
        Mat3 R;                 // body -> global
        Vec3 w, alpha;          // angular velocity and acceleration (global)
    };

    class body_track
    {
    public:
        // bodies: model bodies to record (index k of state() refers to this list)
        void set_bodies(const std::vector<int> &bodies) {bodies_ = bodies; clear();}
        void clear() {t_.clear(); s_.clear();}

        // the state of the model at time t (times increasing)
        void record(const model&, double t);

        int size() const {return (int)t_.size();}
        double t_begin() const {return t_.empty() ? 0.0 : t_.front();}
        double t_end() const {return t_.empty() ? 0.0 : t_.back();}

        // body k at time t (clamped to the recorded interval)
        body_snapshot state(int k, double t) const;

    private:
        std::vector<int> bodies_;
        std::vector<double> t_;
        std::vector<std::vector<body_snapshot>> s_;
    };
}

#endif
