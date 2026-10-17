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
// REEF3D::CORAL: wave kinematics of a resolved fluid (REEF3D FNPF, NHFLOW), sampled at the points of the
// loads that use them (Morison members on bodies and frame members, dynamic mooring lines).
//
// Once per fluid time level the coupling collects the points (load::kinematics_points: strip centres,
// line nodes, at their current positions), samples the fluid there (velocity, free surface elevation,
// whether the point lies in the fluid domain) and adds the sample. Between and after the levels the
// values of a point are a Lagrange polynomial in time through the last order + 1 samples (default 2);
// the local acceleration is its time derivative. A request at an index of an owner returns the values
// of that point; points outside the fluid domain, requests before the first sample and requests without
// an index take the background waves (CORAL's own wave_field, if given) or still water.
// The dynamic pressure is not sampled (0).

#ifndef CORAL_SAMPLED_KINEMATICS_H_
#define CORAL_SAMPLED_KINEMATICS_H_

#include"coral/hydro/coral_waves.h"
#include<deque>
#include<memory>
#include<unordered_map>
#include<vector>

namespace coral
{
    class model;

    class sampled_kinematics : public wave_kinematics
    {
    public:
        explicit sampled_kinematics(std::shared_ptr<const wave_field> background=nullptr) : bg_(std::move(background)) {}

        double elevation(double x, double y, double t) const override;
        void kinematics(const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const override;
        double elevation_at(const void *owner, int index, const Vec3 &p, double t) const override;
        void kinematics_at(const void *owner, int index, const Vec3 &p, double t, Vec3 &u, Vec3 &du, double &p_rho) const override;

        // the points of all loads that take their kinematics from here, at their current positions (in
        // sampling order); a change in the number of points drops the earlier samples
        void collect(const model&, std::vector<Vec3> &points);
        // the fluid at time t at the points of the last collect: velocity, surface elevation, inside (0/1)
        void add_sample(double t, const std::vector<Vec3> &u, const std::vector<double> &eta, const std::vector<int> &inside);

        int npoints() const {return n_;}
        int nsamples() const {return (int)lev_.size();}
        double sample_time() const {return lev_.empty() ? 0.0 : lev_.back().t;}
        const std::shared_ptr<const wave_field> &background() const {return bg_;}

        int order = 2;          // polynomial order in time (0, 1, 2)

    private:
        struct level {double t; std::vector<Vec3> u; std::vector<double> eta; std::vector<int> in;};
        // index of (owner, index) in the samples, -1: none
        int find(const void *owner, int index) const;
        // Lagrange weights of the last levels at t and of their time derivative
        int weights(double t, double *w, double *dw) const;

        std::shared_ptr<const wave_field> bg_;
        std::unordered_map<const void*, std::pair<int,int>> owners_;   // owner -> first point, count
        std::deque<level> lev_;
        int n_ = 0;
    };
}

#endif
