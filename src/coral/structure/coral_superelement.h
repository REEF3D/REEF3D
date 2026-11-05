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
// REEF3D::CORAL: Craig-Bampton superelement of a linear beam substructure (jackets, tripods, monopiles
// in the fixed frame). The members are linear 3D Timoshenko elements (consistent mass, with rotary
// inertia) assembled on nodes merged by position. Member ends on the ground are clamped; ends on rigid
// bodies are interface nodes, rigidly tied to that body (u_node = u_body + theta x r, as SubDyn's rigid
// transition piece). The interior is reduced to the Guyan (static) shapes of the interface motion plus
// n fixed-interface modes:
//   u_I = Phi u_R + Psi q,   Phi = -K_II^-1 K_IR,   K_II psi = w^2 M_II psi (psi^T M_II psi = 1)
//   [M_RR~ M_Rm; M_mR I] [a_R; q''] + [C_RR~ 0; 0 2 zeta Omega] [v_R; q'] + [K_RR~ 0; 0 Omega^2] [u_R; q]
//       = [F_R~; F_m]   (gravity of the substructure, consistent)
// The modal amplitudes q are generalized coordinates of the model (model::add_coordinates).
// Interface motion u_R per body: x - x0 (global) and log(R R0^T) (small rotations, fixed frame);
// forces go to the body as a force at its CoG and a moment. Damping: modal ratio zeta on the modes,
// stiffness-proportional beta on the Guyan part. Analytic tangent (linear in the interface motion).
// The nodes are material points (point_set): recovered positions, velocities and accelerations for
// distributed loads (Morison), forces mapped back through the Guyan and modal shapes.

#ifndef CORAL_SUPERELEMENT_H_
#define CORAL_SUPERELEMENT_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_points.h"
#include"coral/core/coral_types.h"
#include<map>
#include"coral/structure/coral_frame.h"
#include<string>
#include<utility>
#include<vector>

namespace coral
{
    class model;

    // linear 3D Timoshenko beam in local axes (x along the beam): stiffness and consistent mass,
    // order (u, v, w, theta_x, theta_y, theta_z) at end a, then end b
    Eigen::Matrix<double,12,12> beam_stiffness_local(const beam_section&, double L);
    Eigen::Matrix<double,12,12> beam_mass_local(const beam_section&, double L);

    class superelement : public load, public point_set
    {
    public:
        // members: as for the frame builder (body1/body2: -2 free node, -1 ground, >= 0 rigid body);
        // g: gravity acting on the substructure (zero: none); nmode fixed-interface modes are kept
        superelement(model&, const std::string &name, const std::vector<member_spec> &members, int nmode,
                     double zeta, double beta, const Vec3 &g, double tol=1.0e-9);

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {for(int i : bodies_) b.push_back(i);}
        void coupled_coordinates(std::vector<int> &c) const override {for(int k=0; k<nm_; ++k) c.push_back(c0_ + k);}
        void eval(const model&, double t, VecX &f) const override;
        bool has_tangent() const override {return true;}
        void tangent(const model&, double t, double cK, double cC, double cM, std::vector<triplet>&) const override;
        double potential(const model&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        int nmode() const {return nm_;}
        int first_coordinate() const {return c0_;}
        int nnode() const {return (int)p0_.size();}
        int ninterior() const {return nI_;}
        const std::vector<int> &interface_bodies() const {return bodies_;}
        // fixed-interface modal frequencies [Hz]
        const VecX &frequencies() const {return fmode_;}
        // reduced matrices (interface 6 per body, then the modes) and gravity forces
        MatX reduced_mass() const;
        MatX reduced_stiffness() const;
        VecX reduced_gravity() const;
        // current node positions (global, recovered) and the elements as node pairs
        void node_positions(const model&, std::vector<Vec3>&) const;
        const std::vector<std::pair<int,int>> &elements() const {return elem_;}
        // displacement (6 per node: translation, rotation; global) of node k from the current state
        Eigen::Matrix<double,6,1> node_displacement(const model&, int k) const;
        // nodes of a member, in station order (empty if not in this superelement)
        std::vector<int> member_nodes(const std::string &member) const;
        double member_diameter(const std::string &member) const;      // outer diameter of the first section

        // material points: the nodes
        int npoints() const override {return nnode();}
        void point_state(const model&, int k, Vec3 &p, Vec3 &v, Vec3 &a) const override;
        void point_states(const model&, const std::vector<int> &k, std::vector<Vec3> &p, std::vector<Vec3> &v,
                          std::vector<Vec3> &a) const override;
        void add_point_force(const model&, int k, const Vec3 &F, VecX &f) const override;
        void point_coupling(std::vector<int> &bodies, std::vector<int> &coords) const override
        {
            coupled_bodies(bodies);
            coupled_coordinates(coords);
        }

    private:
        struct bstate {Vec3 x, v, a; Mat3 R; Vec3 w, wd;};
        void gather(const model&, std::vector<bstate>&) const;
        // generalized forces: per interface body force (global) and moment (body frame); modes
        void compute(const std::vector<bstate>&, const VecX &q, const VecX &qd, const VecX &qdd, VecX &fb, VecX &fq) const;
        void interface_motion(const std::vector<bstate>&, VecX &u, VecX &v, VecX &a) const;
        void node_state(const model&, const VecX &u, const VecX &vr, const VecX &ar, int k, Vec3 &p, Vec3 &v, Vec3 &a) const;

        std::string name_;
        int nm_ = 0, c0_ = 0, nI_ = 0;
        std::vector<int> bodies_;              // interface bodies
        std::vector<Vec3> xb0_;                // their reference positions
        std::vector<Mat3> Rb0_;                // and orientations
        std::vector<Vec3> p0_;                 // node positions at creation
        std::vector<int> kind_;                // per node: -1 fixed, -2 interior, >= 0 index into bodies_
        std::vector<int> idof_;                // per node: first interior dof (interior nodes)
        std::vector<std::pair<int,int>> elem_;
        std::map<std::string, std::vector<int>> mnodes_;
        std::map<std::string, double> mdiam_;
        MatX KR_, MR_, CR_, MRm_, Phi_, Psi_;  // Guyan stiffness, mass, damping; coupling; recovery
        VecX FR_, Fm_, om2_, cm_, fmode_;      // gravity; Omega^2; modal damping 2 zeta Omega; f [Hz]
    };
}

#endif
