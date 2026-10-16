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
// REEF3D::CORAL: holonomic constraints g(q) = 0 between two bodies (or a body and the ground, index -1),
// and the joints built from them. The Jacobian G (rows of g, columns = generalized velocities of the
// model) gives dg/dt = G nu; the constraint forces on the bodies are -G^T lambda.
//
// Joints are defined in the configuration at their creation: a global point P and a global axis e.
//   fixed        6 equations: point + orientation
//   spherical    3: point
//   revolute     5: point + axis kept parallel           (rotation about e)
//   prismatic    5: orientation + no motion normal to e  (translation along e)
//   cylindrical  4: axis kept parallel + no motion normal to e (translation along and rotation about e)
//   planar       3: no motion along e + axis kept parallel (motion in the plane normal to e)

#ifndef CORAL_CONSTRAINT_H_
#define CORAL_CONSTRAINT_H_

#include"coral/core/coral_types.h"
#include<string>
#include<vector>

namespace coral
{
    class model;

    class constraint
    {
    public:
        virtual ~constraint() = default;
        virtual std::string name() const = 0;
        virtual int size() const = 0;
        // g (size()) and G (size() x ndof), G zero on entry
        virtual void eval(const model&, Eigen::Ref<VecX> g, Eigen::Ref<MatX> G) const = 0;

        int body_a = -1;
        int body_b = -1;
    };

    enum class joint_type {fixed, spherical, revolute, prismatic, cylindrical, planar};

    class joint : public constraint
    {
    public:
        // joint between body a and body b (-1: ground) at the global point P with the global axis e,
        // in the current configuration of the model
        joint(const model&, joint_type, int a, int b, const Vec3 &P, const Vec3 &e, const std::string &name="");

        std::string name() const override {return name_;}
        int size() const override {return (int)prims_.size() + 3*npoint_;}
        void eval(const model&, Eigen::Ref<VecX> g, Eigen::Ref<MatX> G) const override;

        joint_type type() const {return type_;}

        // joint coordinate and its rate: displacement along e (prismatic, cylindrical) or rotation
        // about e (revolute, cylindrical: rotation_coordinate), zero at creation
        double coordinate(const model&) const;
        double rate(const model&) const;
        double rotation_coordinate(const model&) const;
        double rotation_rate(const model&) const;

        // generalized forces of an actuator in the joint: force F along e on b (and -F on a),
        // torque T about e on b (and -T on a)
        void add_axial_force(const model&, double F, VecX &f) const;
        void add_axial_torque(const model&, double T, VecX &f) const;

        // global axis of the joint (fixed to body a)
        Vec3 axis(const model&) const;

    private:
        // primitive conditions
        //   perp: (R_a u) . (R_b w) = 0
        //   dist: (p_b - p_a) . (R_a u) = 0, p = x + R s
        struct prim {int kind; Vec3 u; Vec3 w;};
        void add_perp(const Vec3 &u, const Vec3 &w) {prims_.push_back({0,u,w});}
        void add_dist(const Vec3 &u) {prims_.push_back({1,u,Vec3::Zero()});}

        std::string name_;
        joint_type type_;
        int npoint_ = 0;           // 1: point coincidence p_a = p_b (3 rows)
        Vec3 sa_, sb_;             // joint point in the body frames
        Vec3 ea_, eb_;             // axis in the body frames
        Vec3 pa_, pb_;             // a normal to the axis in the body frames (rotation coordinate)
        double d0_ = 0.0;          // initial distance along the axis
        std::vector<prim> prims_;
    };
}

#endif
