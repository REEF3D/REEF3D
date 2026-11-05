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

#ifndef WAVE_LIB_STEADY_FI_H_
#define WAVE_LIB_STEADY_FI_H_

#include<vector>
#include<functional>

// Velocity potential of a steady periodic wave from its horizontal velocity, for the wave theories that give eta and
// the velocities but no potential (cnoidal, shallow-water, Stokes 5th of Skjelbreia & Hendrickson). FNPF relaxes
// Fifsf to the potential in the generation zone (and starts from it, B 115 / fullini); without it the zone relaxed
// Fifsf to an undefined value and absorbed most of the wave it generated.
//
//   fi(x,z,t) = F(x - c t, z) + b t
//
// F(s,z) is the zero-mean antiderivative along the direction of travel of the horizontal velocity u(s,z) over one
// wave length L (the mean current is left out), computed spectrally from N samples per wave length at nz levels between
// the bed and the highest crest and interpolated (periodic cubic in s, cubic in z). b is the Bernoulli constant at the
// surface, b = < c u - (u^2 + w^2)/2 - g eta > over one wave length, so that fi_t + |grad fi|^2/2 + g eta = 0 holds at
// the surface on average (the FNPF free-surface condition); without it Fifsf of FNPF and the target drift apart in time.
// The table is built once, at the first call, from the wave at that time.

class wave_lib_steady_fi
{
public:
    // L wave length, c phase speed (waves travel in +s), h water depth, g gravity, t0 time of the functions below,
    // uh(s,z), w(s,z), eta(s): the wave at time t0 along the direction of travel (s = x of the wave frame)
    void build(double L, double c, double h, double g, double t0,
               const std::function<double(double,double)> &uh,
               const std::function<double(double,double)> &w,
               const std::function<double(double)> &eta);

    bool built() const {return ok;}

    double fi(double s, double z, double t) const;

private:
    bool ok=false;
    int N=256, nz=33;
    double L=1.0, c=0.0, t0=0.0, b=0.0, zb=0.0, zt=0.0, dz=1.0;
    std::vector<double> F;   // F[kz*N + i] at s_i = i L/N, z_kz = zb + kz dz

    double row(int kz, double s) const;
};

#endif
