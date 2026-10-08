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
// REEF3D::CFD and REEF3D::NHFLOW.   ctrl.txt: Z 30 1, Z 31 dt_print;  input file fem.dat
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
//
// NHFLOW (fem_coupling_nhflow.cpp): the same algorithm on the sigma grid.
// Velocities are sampled and the forcing is spread with the kernel over the
// horizontal cell centres and, vertically, over the cell centres of the wet
// column (renormalised at the bed and the free surface); the forcing is added
// to U,V,W and UH,VH,WH. The probed pressure is the non-hydrostatic pressure
// (filtered over a few steps) plus the hydrostatic pressure below the local free
// surface. The cells inside deformable structures are marked solid (p->DF -1,
// d->solid_flux): no fluxes through the structure; a structure standing dry
// closes its columns over the full height (no overtopping in NHFLOW).

#include"increment.h"
#include"fem_solid.h"
#include<vector>
#include<string>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class field;
class slice;

class fem_coupling : public increment
{
public:
    fem_coupling(lexer*, ghostcell*);
    virtual ~fem_coupling();

    // CFD: adds the forcing (acceleration) to the staggered fields fx, fy, fz
    void start_cfd(lexer*, fdm*, ghostcell*, double alpha, field&, field&, field&, field&, field&, field&, bool finalize);

    // NHFLOW: adds the forcing to U,V,W and UH,VH,WH of the stage and marks the
    // cells inside the structure as solid (p->DF -1: no fluxes through its faces).
    // solid: immersed solids of the grid in d->SOLID; dfreset: nhflow_forcing has
    // reset p->DF in this stage. Ghost cells of the velocities are updated by nhflow_forcing.
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, double alpha, double *UH, double *VH, double *WH, slice &WL, bool solid, bool dfreset, bool finalize);

private:
    struct lpoint
    {
        int face;
        double s, t;        // bilinear position on the face
        double frac;        // fraction of the face area
    };

    static constexpr size_t BP = 23;     // buffer entries per Lagrangian point / debris particle

    void ini_points(lexer*, ghostcell*);
    void point_state(int q, fem_solid::Vec3& xp, fem_solid::Vec3& vp, fem_solid::Vec3& n, double& A) const;
    double interpolate_kernel(lexer*, field&, double, double, double, int comp, field *ro=nullptr, double *rho=nullptr);
    void spread(lexer*, field&, field&, field&, const fem_solid::Vec3& xp, const fem_solid::Vec3& f, double A, const fem_solid::Vec3* n);
    double kernel(double) const;
    void finish_step(lexer*, ghostcell*, double alpha);
    void sample_bed(lexer*, fdm*, ghostcell*);
    void probe_pressure(lexer*, fdm*, const fem_solid::Vec3& xp, const fem_solid::Vec3& n, double *b, bool hydrostatic=false);
    bool probe_fallback(int q, const fem_solid::Vec3& xp, fem_solid::Vec3& fb) const;   // probe beside a rigid body, at the height of the point
    void probe_beside(lexer*, fdm*, const fem_solid::Vec3& xp, fem_solid::Vec3 q, double *b);
    std::vector<fem_solid::Vec3> rb_lo, rb_hi;                // horizontal bounding boxes of the rigid bodies
    void rigid_boxes();
    bool on_bed(int q, const fem_solid::Vec3& nf, double zmin) const;   // rigid face near the bed: no forcing of the water under it
    void pressure_loads(lexer*, ghostcell*, std::vector<fem_solid::Vec3>& F, bool hydrostatic=false);
    void print(lexer*);
    void first_call(lexer*, ghostcell*);     // supports on the bed, check mode, gravity settling

    // fluid of the current call: CFD (cfd) or NHFLOW (nhf, WLn, nhf_solid)
    fdm *cfd = nullptr;
    fdm_nhf *nhf = nullptr;
    slice *WLn = nullptr;
    bool nhf_solid = false;
    bool nhflow = false;

    // NHFLOW (fem_coupling_nhflow.cpp)
    double nhf_column(lexer*, double *f, int i, int j, double z, bool faces) const;
    double nhf_ipol(lexer*, double *f, double x, double y, double z, bool faces) const;
    double nhf_surface(lexer*, double x, double y) const;
    double nhf_bedlevel(lexer*, double x, double y) const;
    bool nhf_wet_column(lexer*, double x, double y) const;
    double nhf_level(lexer*, const fem_solid::Vec3& x) const;  // signed distance to the bed and the solids of the grid
    bool nhf_in_solid(lexer*, const fem_solid::Vec3& x) const;  // bed, solids of the grid, cells marked solid
    double nhf_hv(lexer*, int i, int j, double z) const;       // vertical kernel width at the point
    double nhf_kernel_ipol(lexer*, double *f, const fem_solid::Vec3& x) const;
    void nhf_spread(lexer*, const fem_solid::Vec3& xp, const fem_solid::Vec3& du, double A, const fem_solid::Vec3* n);
    void nhf_apply_forcing(lexer*, double *UH, double *VH, double *WH);   // summed increments, at most one full correction per cell
    std::vector<double> sp_w, sp_u;
    void nhf_probe_pressure(lexer*, const fem_solid::Vec3& xp, const fem_solid::Vec3& n, double *b);
    void nhf_probe_beside(lexer*, const fem_solid::Vec3& xp, fem_solid::Vec3 q, double *b);
    void nhf_sample_bed(lexer*, ghostcell*);
    void nhf_mark_solid(lexer*, ghostcell*, bool dfreset);    // cells inside the deformable structure: p->DF = -1
    std::vector<int> df_cells, df_save;
    bool warned_overtop = false;
    // factor on the added-mass estimate of rigid body k: given, or 1 (CFD); NHFLOW: 5 (the
    // non-hydrostatic pressure responds more strongly to the forcing of a body; with 1, light
    // debris in a shallow, fast flow became unstable), but at most 30 body masses and at least
    // 2 (a much larger added mass delays the response of very light bodies such as empty
    // containers: the impulse of a slam was stored and given back, the body was thrown out;
    // see also nhf_limit_rigid)
    double added_mass_factor(int k) const;
    std::vector<double> pnh_bar, pnh_bar_b;    // NHFLOW: filtered non-hydrostatic pressure at the probes
    void nhf_filter_pressure();
    void nhf_limit_rigid(lexer*, bool apply);     // rigid bodies: speed relative to the water at most the terminal speed
    std::vector<fem_solid::Vec3> lim_u;
    std::vector<double> lim_vt;
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
        double sutil = -1.0, t_sutil = 0.0;          // reinforcement
        double t_syield = -1.0;
        fem_solid::Vec3 x_syield = fem_solid::Vec3::Zero();
        int nrupt = 0;
    } sm;

    double printtime;
    int printcount;
    double starttime;
    std::string outdir;
};

#endif
