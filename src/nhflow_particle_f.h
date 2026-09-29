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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_PARTICLE_F_H_
#define NHFLOW_PARTICLE_F_H_

#include"nhflow_particle.h"
#include"increment.h"
#include<vector>
#include<random>
#include<fstream>
#include<cstdint>

class slice;

// particle states, stored as double in nhflow_particle_data::state
constexpr double NHFP_WATER    = 0.0;  // in the water column
constexpr double NHFP_SURFACE  = 1.0;  // at the free surface (floating, subject to windage)
constexpr double NHFP_STRANDED = 2.0;  // in a dry cell
constexpr double NHFP_BED      = 3.0;  // deposited on the bed

struct nhflow_particle_data
{
    double x,y,z;          // position, model coordinates
    double u0,v0,w0;       // drift velocity at (x^n,t^n) for the Heun step
    double u,v,w;          // fluid velocity at the particle, output
    double ws;             // rise (>0) or settling (<0) velocity [m/s]
    double cw;             // windage coefficient [-], fraction of U10
    double mode;           // 0: 3D particle, 1: surface-trapped
    double state;          // nhflow_particle_state
    double id;             // global id, independent of the decomposition
    double src;            // release number
    double t0;             // release time
};

class nhflow_particle_f : public nhflow_particle, public increment
{
public:
    nhflow_particle_f(lexer*, ghostcell*);
    virtual ~nhflow_particle_f();

    void ini(lexer*, fdm_nhf*, ghostcell*) override;
    void step_begin(lexer*, fdm_nhf*, ghostcell*) override;
    void step_end(lexer*, fdm_nhf*, ghostcell*) override;

private:
    struct release
    {
        int type;            // 1 point, 2 line, 3 box
        double xs,ys,zs,xe,ye,ze;
        long long num;
        double ts,te;
        double ws,cw;
        int mode;
        long long idoffset;
        long long released;
    };

    // seeding
    void setup_releases(lexer*);
    void seed(lexer*, fdm_nhf*, ghostcell*, double, bool);
    double release_time(const release&, long long) const;
    void release_position(const release&, long long, double&, double&, double&) const;
    static double hash_uniform(uint64_t, uint64_t);

    // motion
    void drift_velocity(lexer*, fdm_nhf*, nhflow_particle_data&, double&, double&, double&);
    void fluid_velocity(lexer*, fdm_nhf*, double, double, double, double&, double&, double&);
    void column_velocity(lexer*, fdm_nhf*, int, int, double, double&, double&, double&);
    double column_scalar(lexer*, double*, int, int, double);
    double scalar_ipol(lexer*, double*, double, double, double);
    double diffusivity(lexer*, fdm_nhf*, double, double, double);
    void vertical_diffusion(lexer*, fdm_nhf*, nhflow_particle_data&);
    void domain_boundary(lexer*);
    void vertical_bounds(lexer*, fdm_nhf*, nhflow_particle_data&);
    void wetdry(lexer*, fdm_nhf*, nhflow_particle_data&);
    double surface(lexer*, fdm_nhf*, double, double);
    double bedlevel(lexer*, fdm_nhf*, double, double);
    bool owned(lexer*, double, double) const;

    // parallel
    void xchange(lexer*, ghostcell*);

    // output
    void print(lexer*, fdm_nhf*, ghostcell*);
    void gather(lexer*, ghostcell*, std::vector<nhflow_particle_data>&);
    void print_vtp(lexer*, const std::vector<nhflow_particle_data>&);
    void print_csv(lexer*, const std::vector<nhflow_particle_data>&);
    void print_log(lexer*, ghostcell*);

    std::vector<nhflow_particle_data> P;
    std::vector<release> R;

    std::mt19937_64 rng;
    std::normal_distribution<double> gauss;

    double xloc_s,xloc_e,yloc_s,yloc_e;
    double U10,cosw,sinw;
    double printtime;
    int printcount;
    long long numout,numout_global;
    std::ofstream logout;
    std::ofstream csvout;

    static const int NF = sizeof(nhflow_particle_data)/sizeof(double);
};

#endif
