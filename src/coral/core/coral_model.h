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
// REEF3D::CORAL: the multibody model: bodies, generalized coordinates, constraints and loads.
// Generalized velocities: 6 per body, [v (global); w (body frame)], body i at 6*i, followed by the
// generalized coordinates (scalar degrees of freedom such as the modal amplitudes of a superelement),
// coordinate k at 6*nbody() + k. Their inertia comes from the loads that own them.
// Body index -1 is the ground (fixed, x = 0, R = I).

#ifndef CORAL_MODEL_H_
#define CORAL_MODEL_H_

#include"coral/core/coral_types.h"
#include"coral/core/coral_body.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_load.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    class model
    {
    public:
        model() = default;
        model(const model&) = delete;
        model &operator=(const model&) = delete;

        // the body's current configuration becomes its reference configuration (x0, R0)
        int add_body(const rigid_body&);
        int find_body(const std::string&) const;     // -1 if not found
        // mass at the global point p rigidly attached to body i (inertia J about its own centre, global
        // axes): part of the body's inertia forces and of its weight (gravity_load)
        void add_point_mass(int i, const Vec3 &p, double mass, const Mat3 &J=Mat3::Zero());
        constraint &add_constraint(std::unique_ptr<constraint>);
        load &add_load(std::unique_ptr<load>);

        // n generalized coordinates (value, rate and acceleration zero); returns the index of the first
        int add_coordinates(int n, const std::string &name);
        int ncoord() const {return (int)cq.size();}
        int coord_dof(int k) const {return 6*nbody() + k;}

        int nbody() const {return (int)bodies.size();}
        int ndof() const {return 6*nbody() + ncoord();}
        int ncon() const;

        // body i, or the ground for i = -1
        const rigid_body &body(int i) const {return i<0 ? ground_ : bodies[i];}

        // loads, inertia forces (M dnu/dt + gyroscopic terms) and constraints of the current state;
        // fd_only: only the loads without an analytic tangent (finite-difference Jacobians)
        void forces(double t, VecX &f, bool fd_only=false) const;
        void inertia_forces(VecX &r) const;
        void constraints(VecX &g, MatX &G) const;

        void initialize(double t0);
        void commit(double t, double dt);

        // helpers for loads: force F (global) at the global point p of body i, moment M (global)
        void add_force(VecX &f, int i, const Vec3 &p, const Vec3 &F) const;
        void add_moment(VecX &f, int i, const Vec3 &M) const;

        // analytic tangents of the loads that have one: cK dF/dq + cC dF/dnu + cM dF/dnu'
        bool has_tangents() const;
        void tangents(double t, double cK, double cC, double cM, std::vector<triplet> &T) const;

        // distance-2 colouring of the coupling graph (loads and constraints) of the nodes: bodies
        // 0..nbody()-1, then the generalized coordinates. Nodes of the same colour share no load or
        // constraint and no neighbour, so a finite-difference Jacobian can perturb them together.
        // Returns the number of colours; nbr[i]: node i and its neighbours (sorted). fd_only: the
        // loads with an analytic tangent are left out.
        int coloring(std::vector<int> &color, std::vector<std::vector<int>> &nbr, bool fd_only=false) const;
        // degrees of freedom of node i of the colouring (6 for a body, 1 for a coordinate)
        int node_dofs(int i) const {return i < nbody() ? 6 : 1;}
        int node_dof0(int i) const {return i < nbody() ? 6*i : 6*nbody() + (i - nbody());}

        // configuration (positions, rotations, coordinates) and increments of it:
        // x = x0 + s dx, R = R0 exp(s dtheta), cq = cq0 + s dq
        struct config {std::vector<Vec3> x; std::vector<Mat3> R; VecX cq;};
        config save_config() const;
        void set_config(const config&, const VecX &dq, double s=1.0);
        // velocities and accelerations of bodies and coordinates to zero
        void zero_rates();
        // largest entry of an increment of the generalized coordinates (6 per body, then coordinates),
        // the coordinates weighted with cscale
        double increment_norm(const VecX &d) const;

        double kinetic_energy() const;
        // restart from a checkpoint: the state of the bodies (and their mass and inertia, which the
        // fluid solver may set), of the generalized coordinates and of the loads with a history
        void save_state(std::vector<double>&) const;
        bool load_state(const std::vector<double>&, size_t &k);
        double potential_energy() const;

        Vec3 gravity = Vec3(0.0, 0.0, -9.81);

        std::vector<rigid_body> bodies;
        std::vector<std::unique_ptr<constraint>> constraint_list;
        std::vector<std::unique_ptr<load>> load_list;

        struct point_mass {int body; Vec3 s; double mass; Mat3 J;};    // s, J in the body frame
        std::vector<point_mass> point_masses;

        // generalized coordinates: value, rate, acceleration; displacement per unit coordinate (for the
        // convergence tests of the solvers, set by the owning load, default 1); names of the blocks
        VecX cq, cv, ca, cscale;
        std::vector<std::pair<std::string,int>> coord_blocks;     // name, first index

    private:
        rigid_body ground_;
    };
}

#endif
