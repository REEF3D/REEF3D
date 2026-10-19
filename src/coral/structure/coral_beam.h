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
// REEF3D::CORAL: co-rotational beam element between two attachments (a node body, a rigid body at an
// offset, or the ground). Large displacements and rotations, small strains: the element frame
// follows the chord and the mean of the end triads (Crisfield, Battini); in it the deformations are
// the elongation u = L - L0, the torsion and the bending rotations of both ends against the chord,
// with the linear Timoshenko stiffness of the section:
//   N = EA u/L0,  T = GJ (tb_x - ta_x)/L0,
//   M_a = EI/(L0 (1+P)) ((4+P) ta + (2-P) tb),  M_b = EI/(L0 (1+P)) ((2-P) ta + (4+P) tb),
//   P = 12 EI/(GA L0^2) (0 for GA = 0: Euler-Bernoulli), about the local y and z axes.
// The generalized forces on the attached bodies are -B^T f with B = dd/dq analytic (q: position global,
// rotation body frame of each attached body): the spin of the element frame follows from the chord and
// the mean y axis of the end triads, the end rotations through the inverse left Jacobian of SO(3).
// Stiffness-proportional damping beta K. Analytic tangent: material and geometric stiffness by
// central differences of -B^T f (with the analytic B), damping -beta B^T K B.
// The mass is lumped on the nodes by the frame builder (coral_frame).
// Local axes: x along the chord, y from the orientation vector at creation, z = x cross y.

#ifndef CORAL_BEAM_H_
#define CORAL_BEAM_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_types.h"

namespace coral
{
    struct beam_section
    {
        double EA = 0.0, EIy = 0.0, EIz = 0.0, GJ = 0.0;
        double GAy = 0.0, GAz = 0.0;             // shear stiffness (0: no shear deformation)
        double m = 0.0;                          // mass per length [kg/m]
        double rIx = 0.0, rIy = 0.0, rIz = 0.0;  // rotary inertia per length (rho J, rho Iy, rho Iz) [kg m]
        double D = 0.0;                          // outer diameter (hydrodynamics, output)

        // circular tube, outer diameter D, wall t; shear factor 0.5 of a thin tube
        static beam_section tube(double D, double t, double E, double G, double rho);
        static beam_section blend(const beam_section &a, const beam_section &b, double s);    // (1-s) a + s b
    };

    struct beam_end
    {
        int body = -1;               // -1: ground
        Vec3 p = Vec3::Zero();       // global position at creation
    };

    class beam_element : public load
    {
    public:
        beam_element(const model&, const beam_end &a, const beam_end &b, const beam_section&,
                     const Vec3 &orient=Vec3::UnitZ(), double beta=0.0, const std::string &name="beam");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(ba_); b.push_back(bb_);}
        void eval(const model&, double t, VecX &f) const override;
        bool has_tangent() const override {return true;}
        void tangent(const model&, double t, double cK, double cC, double cM, std::vector<triplet>&) const override;
        double potential(const model&) const override;

        // B = dd/dq of the current state (7 x 12: dx_a, dtheta_a, dx_b, dtheta_b), for tests
        Eigen::Matrix<double,7,12> strain_matrix(const model&) const;

        double length0() const {return L0_;}
        const beam_section &section() const {return s_;}
        // deformations d = (u, ta_x, ta_y, ta_z, tb_x, tb_y, tb_z) and local forces f = K d
        Eigen::Matrix<double,7,1> deformation(const model&) const;
        Eigen::Matrix<double,7,1> local_forces(const Eigen::Matrix<double,7,1> &d) const;
        // current element frame (columns x, y, z) and end positions
        Mat3 frame(const model&, Vec3 &pa, Vec3 &pb) const;

    private:
        using Vec7 = Eigen::Matrix<double,7,1>;
        using Mat7x12 = Eigen::Matrix<double,7,12>;
        using Vec12 = Eigen::Matrix<double,12,1>;
        Eigen::Matrix<double,7,1> deformation(const Vec3 &pa, const Mat3 &Ta, const Vec3 &pb, const Mat3 &Tb) const;
        void ends(const model&, Vec3 &pa, Mat3 &Ta, Vec3 &pb, Mat3 &Tb) const;
        // deformation and B for end bodies at (xa, Ra), (xb, Rb)
        void kinematics(const Vec3 &xa, const Mat3 &Ra, const Vec3 &xb, const Mat3 &Rb, Vec7 &d, Mat7x12 &B) const;
        // element forces on the 12 coordinates for the given end states (velocities v global, w body)
        Vec12 element_force(const Vec3 &xa, const Mat3 &Ra, const Vec3 &xb, const Mat3 &Rb, const Vec12 &qd) const;
        Eigen::Matrix<double,7,7> local_stiffness() const;
        Vec12 rates(const model&) const;

        std::string name_;
        int ba_, bb_;
        Vec3 sa_, sb_;           // end points in the body frames
        Mat3 Qa_, Qb_;           // element frame at creation in the body frames
        double L0_;
        beam_section s_;
        double beta_;
    };
}

#endif
