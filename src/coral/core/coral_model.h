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
// REEF3D::CORAL: the multibody model: bodies, constraints and loads.
// Generalized velocities: 6 per body, [v (global); w (body frame)], body i at 6*i.
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
        constraint &add_constraint(std::unique_ptr<constraint>);
        load &add_load(std::unique_ptr<load>);

        int nbody() const {return (int)bodies.size();}
        int ndof() const {return 6*nbody();}
        int ncon() const;

        // body i, or the ground for i = -1
        const rigid_body &body(int i) const {return i<0 ? ground_ : bodies[i];}

        // loads, inertia forces (M dnu/dt + gyroscopic terms) and constraints of the current state
        void forces(double t, VecX &f) const;
        void inertia_forces(VecX &r) const;
        void constraints(VecX &g, MatX &G) const;

        void initialize(double t0);
        void commit(double t, double dt);

        // helpers for loads: force F (global) at the global point p of body i, moment M (global)
        void add_force(VecX &f, int i, const Vec3 &p, const Vec3 &F) const;
        void add_moment(VecX &f, int i, const Vec3 &M) const;

        // distance-2 colouring of the coupling graph of the bodies (loads and constraints): bodies of the
        // same colour share no load or constraint and no neighbour, so a finite-difference Jacobian can
        // perturb them together. Returns the number of colours; nbr[i]: body i and its neighbours.
        int coloring(std::vector<int> &color, std::vector<std::vector<int>> &nbr) const;

        double kinetic_energy() const;
        double potential_energy() const;

        Vec3 gravity = Vec3(0.0, 0.0, -9.81);

        std::vector<rigid_body> bodies;
        std::vector<std::unique_ptr<constraint>> constraint_list;
        std::vector<std::unique_ptr<load>> load_list;

    private:
        rigid_body ground_;
    };
}

#endif
