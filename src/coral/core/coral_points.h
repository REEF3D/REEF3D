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
// REEF3D::CORAL: material points of the structure whose motion follows the generalized coordinates, for
// distributed loads that need positions, velocities and accelerations along a structure and give back
// point forces (Morison loads on members, later aerodynamic and ice loads). Implemented by rigid-body
// points (frame nodes, points on bodies) and by superelements (recovered node motion).

#ifndef CORAL_POINTS_H_
#define CORAL_POINTS_H_

#include"coral/core/coral_types.h"
#include<vector>

namespace coral
{
    class model;

    class point_set
    {
    public:
        virtual ~point_set() = default;
        virtual int npoints() const = 0;
        // position, velocity, acceleration (global) of point k in the current state
        virtual void point_state(const model&, int k, Vec3 &p, Vec3 &v, Vec3 &a) const = 0;
        // the same for several points at once (implementations may share work between the points)
        virtual void point_states(const model &m, const std::vector<int> &k, std::vector<Vec3> &p, std::vector<Vec3> &v,
                                  std::vector<Vec3> &a) const
        {
            p.resize(k.size()); v.resize(k.size()); a.resize(k.size());
            for(size_t i=0; i<k.size(); ++i) point_state(m, k[i], p[i], v[i], a[i]);
        }
        // generalized forces of the force F (global) acting at point k
        virtual void add_point_force(const model&, int k, const Vec3 &F, VecX &f) const = 0;
        // bodies and generalized coordinates the points move with
        virtual void point_coupling(std::vector<int> &bodies, std::vector<int> &coords) const = 0;
    };

    // points fixed to rigid bodies (frame nodes: the body's centre)
    class body_points : public point_set
    {
    public:
        // point at the global position p on body b (current configuration)
        int add(const model&, int b, const Vec3 &p);
        int npoints() const override {return (int)body_.size();}
        void point_state(const model&, int k, Vec3 &p, Vec3 &v, Vec3 &a) const override;
        void add_point_force(const model&, int k, const Vec3 &F, VecX &f) const override;
        void point_coupling(std::vector<int> &bodies, std::vector<int>&) const override {for(int b : body_) bodies.push_back(b);}

    private:
        std::vector<int> body_;
        std::vector<Vec3> s_;
    };
}

#endif
