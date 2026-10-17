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
// REEF3D::CORAL inside REEF3D: sampling of the resolved fluid at points for CORAL's sampled wave
// kinematics (fluid_kinematics in coral.txt; coral/hydro/coral_sampled_kinematics.h).
// Each rank fills the points that lie in its own cells (columns of level 0), the caller sums over the
// ranks. Per point (REEF3D frame): u, v, w, z of the free surface, and 1 if the point was sampled
// (0 elsewhere). Above the free surface the velocity of the surface is taken, below the bed zero.

#ifndef CORAL_FLUID_SAMPLER_H_
#define CORAL_FLUID_SAMPLER_H_

#include<vector>

class lexer;
class fdm_fnpf;
class fdm_nhf;

class coral_fluid_sampler
{
public:
    virtual ~coral_fluid_sampler() = default;
    // xyz: 3 per point; out: 5 per point (u, v, w, surface z, sampled), filled for the own points only
    virtual void sample(const std::vector<double> &xyz, std::vector<double> &out) = 0;
    virtual const char* name() const = 0;
};

// FNPF: velocities of the sigma grid (U, V, W of the last velcalc), free surface WL + bed
class coral_sampler_fnpf : public coral_fluid_sampler
{
public:
    coral_sampler_fnpf(lexer *p, fdm_fnpf *c) : p_(p), c_(c) {}
    void sample(const std::vector<double>&, std::vector<double>&) override;
    const char* name() const override {return "FNPF";}
private:
    lexer *p_;
    fdm_fnpf *c_;
};

// NHFLOW: cell-centred velocities of the sigma grid, free surface WL + bed
class coral_sampler_nhflow : public coral_fluid_sampler
{
public:
    coral_sampler_nhflow(lexer *p, fdm_nhf *d) : p_(p), d_(d) {}
    void sample(const std::vector<double>&, std::vector<double>&) override;
    const char* name() const override {return "NHFLOW";}
private:
    lexer *p_;
    fdm_nhf *d_;
};

#endif
