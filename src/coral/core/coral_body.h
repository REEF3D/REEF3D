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
// REEF3D::CORAL: rigid body.
// Configuration: x (centre of gravity, global) and R (body -> global).
// Generalized velocity nu = [v; w]: v velocity of the CoG in the global frame, w angular velocity in
// the body frame. Generalized force = [F (global); M about the CoG (body frame)].
// Points on the body are given by s, the position relative to the CoG in the body frame.

#ifndef CORAL_BODY_H_
#define CORAL_BODY_H_

#include"coral/core/coral_types.h"
#include<string>

namespace coral
{
    struct rigid_body
    {
        std::string name;
        double mass = 0.0;
        Mat3 inertia = Mat3::Zero();    // about the CoG, body frame

        Vec3 x = Vec3::Zero();
        Mat3 R = Mat3::Identity();
        Vec3 v = Vec3::Zero();
        Vec3 w = Vec3::Zero();
        Vec3 vd = Vec3::Zero();         // dv/dt (global)
        Vec3 wd = Vec3::Zero();         // dw/dt (body frame)

        // configuration when the body was added to the model (reference of linear models and joints)
        Vec3 x0 = Vec3::Zero();
        Mat3 R0 = Mat3::Identity();

        Vec3 point(const Vec3 &s) const {return x + R*s;}
        Vec3 point_velocity(const Vec3 &s) const {return v + R*w.cross(s);}
        Vec3 point_acceleration(const Vec3 &s) const {return vd + R*(wd.cross(s) + w.cross(w.cross(s)));}
        Vec3 to_body(const Vec3 &p) const {return R.transpose()*(p - x);}
        Vec3 omega_global() const {return R*w;}
        Vec3 angular_momentum() const {return R*(inertia*w);}    // global, about the CoG

        double kinetic_energy() const {return 0.5*mass*v.squaredNorm() + 0.5*w.dot(inertia*w);}
    };
}

#endif
