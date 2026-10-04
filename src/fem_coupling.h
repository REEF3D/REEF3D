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

#ifndef FEM_COUPLING_H_
#define FEM_COUPLING_H_

// Two-way coupling of the explicit FEM solid solver (fem_solid) with
// REEF3D::CFD.   ctrl.txt: Z 30 1, Z 31 dt_print;  input file fem.dat
//
// Per RK stage:
//   - Lagrangian points on the surface faces of the intact elements carry
//     the structural velocity; direct forcing (u_s - u)/(alpha dt) is spread
//     with the Roma kernel onto the staggered forcing fields (as FSI strips)
//   - debris particles (nodes without intact elements) get a quadratic drag,
//     the reaction is spread back onto the fluid
// Final stage:
//   - loads reaction: the fluid parcel in the forcing volume of every
//     surface point (rho dV, u_f) is attached to the face nodes for the step
//     and moves with them through the substeps, its momentum exchange is the
//     load; the enclosed fluid gives buoyancy and inertia per node (m - rho_f V)
//   - loads pressure: the pressure is probed outside every surface point and
//     integrated as -p n dA onto the face nodes (explicit)
//   - loads hybrid (default): reaction, plus a low-pass filtered correction
//     (pressure - reaction) per node: stable like the parcels, mean loads from
//     the probed pressure (the fluid inside the immersed boundary distorts the
//     parcel loads of one-sided and partly wetted structures)
//   - debris: drag + buoyancy
//   - the solid is advanced over the fluid time step with subcycling
//
// The solid is replicated on every rank (owner ranks sample, one
// MPI_Allreduce per stage), all ranks advance it identically.

#include"increment.h"
#include"fem_solid.h"
#include<vector>
#include<string>

class lexer;
class fdm;
class ghostcell;
class field;

class fem_coupling : public increment
{
public:
    fem_coupling(lexer*, ghostcell*);
    virtual ~fem_coupling();

    // CFD: adds the forcing (acceleration) to the staggered fields fx, fy, fz
    void start_cfd(lexer*, fdm*, ghostcell*, double alpha, field&, field&, field&, field&, field&, field&, bool finalize);

private:
    struct lpoint
    {
        int face;
        double s, t;        // bilinear position on the face
        double frac;        // fraction of the face area
    };

    void ini_points(lexer*, ghostcell*);
    void point_state(int q, fem_solid::Vec3& xp, fem_solid::Vec3& vp, fem_solid::Vec3& n, double& A) const;
    double interpolate_kernel(lexer*, field&, double, double, double, int comp, field *ro=nullptr, double *rho=nullptr);
    void spread(lexer*, field&, field&, field&, const fem_solid::Vec3& xp, const fem_solid::Vec3& f, double A, const fem_solid::Vec3* n);
    double kernel(double) const;
    void finish_step(lexer*, ghostcell*, double alpha);
    void probe_pressure(lexer*, fdm*, const fem_solid::Vec3& xp, const fem_solid::Vec3& n, double *b, bool hydrostatic=false);
    void pressure_loads(lexer*, fdm*, ghostcell*, std::vector<fem_solid::Vec3>& F, bool hydrostatic=false);
    void print(lexer*);
    void first_call(lexer*, fdm*, ghostcell*);     // supports on the bed, check mode, gravity settling
    void warnings(lexer*, std::vector<std::string>&);
    void update_summary(lexer*, bool write);
    void write_summary(lexer*);

    fem_solid fs;

    std::vector<lpoint> pts;
    int surf_version;
    double dxmin;

    std::vector<double> buf;        // sampled fluid data, see start_cfd
    std::vector<fem_solid::Vec3> fdeb;   // debris drag of the current stage
    std::vector<fem_solid::Vec3> fprb, fp_bar, fr_bar;   // hybrid loads: probed pressure, filtered pressure and parcel loads
    bool hybrid_ini = false;

    double rho_w;
    bool initialised;
    double force_scale;             // 1, or 1/slice width in 2D (forces per metre)
    int nstep;

    // engineering summary (running maxima and events)
    struct summary
    {
        double shear = 0.0, t_shear = 0.0;
        double moment = 0.0, t_moment = 0.0;
        double fluid = 0.0, t_fluid = 0.0;
        double disp = 0.0, t_disp = 0.0;
        double util = -1.0, t_util = 0.0;
        fem_solid::Vec3 x_util = fem_solid::Vec3::Zero();
        double t_crack = -1.0;
        fem_solid::Vec3 x_crack = fem_solid::Vec3::Zero();
        double t_fail = -1.0;
        fem_solid::Vec3 x_fail = fem_solid::Vec3::Zero();
        bool yielding = false;
    } sm;

    double printtime;
    int printcount;
    double starttime;
    std::string outdir;
};

#endif
