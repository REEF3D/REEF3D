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

/*--------------------------------------------------------------------
dem_f: REEF3D adapter of the DEM core.

Coupling modes (ctrl E 11):
 0 dry        no fluid forces
 1 unresolved particles smaller than the grid: drag (Haider-Levenspiel with
              Di Felice voidage correction, implicit in the DEM step), buoyancy from
              volume quadrature, added mass and fluid acceleration; the reaction is
              spread to the fluid momentum with a compact kernel
 2 resolved   particles larger than the grid: direct forcing of the fluid velocity
              inside the particle level set, hydrodynamic force from the volume
              integral of the forcing (Uhlmann 2005)
 3 hybrid     per particle, resolved if d_eq/dx >= E 12, else unresolved

The particle state is replicated on all ranks. Fluid data and forces are
evaluated on the rank that owns the evaluation point and summed with one
MPI_Allreduce per step. Rank 0 broadcasts the state after every step.
--------------------------------------------------------------------*/

#ifndef DEM_F_H_
#define DEM_F_H_

#include"dem.h"
#include"dem_core.h"
#include"increment.h"
#include<vector>
#include<string>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class field;
class field1;
class field2;
class field3;
class field4a;
class slice;

using namespace std;

class dem_f final : public dem, public increment
{
public:
    dem_f(lexer*, ghostcell*);
    virtual ~dem_f();

    void start_cfd(lexer*, fdm*, ghostcell*) override final;
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*) override final;

    void forcing_cfd(lexer*, fdm*, ghostcell*, int, double, field&, field&, field&, bool) override final;
    void forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, int, double, double*, double*, double*, slice&, bool, bool) override final;

private:
    // setup
    void read_input(lexer*, ghostcell*);
    void ini(lexer*, ghostcell*);
    void ini_cfd(lexer*, ghostcell*);
    void ini_nhflow(lexer*, ghostcell*);

    // time step
    void step(lexer*, ghostcell*);
    void set_loads(lexer*, ghostcell*);
    void sync(lexer*, ghostcell*);
    void deactivate(lexer*);
    int substeps(lexer*);

    // point ownership in the domain decomposition
    bool owns(lexer*, double, double, double);

    // CFD coupling
    void fluid_cfd(lexer*, fdm*, ghostcell*);
    void feedback_cfd(lexer*, fdm*, ghostcell*);
    void walls_cfd(lexer*, fdm*, ghostcell*, double, vector<dem_contact>&);
    double wallphi_cfd(lexer*, fdm*, double, double, double);

    // NHFLOW coupling
    void fluid_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void feedback_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void walls_nhflow(lexer*, fdm_nhf*, ghostcell*, double, vector<dem_contact>&);
    double wallphi_nhflow(lexer*, fdm_nhf*, double, double, double);

    // common
    void gather_walls(lexer*, ghostcell*, const vector<double>&, vector<dem_contact>&);
    void drag(lexer*, int, double, double, double);
    double heaviside(double, double);

    // grid point relative to a particle centre; in 2D the y offset is ignored (x-z slab)
    dem_vec relpos(lexer*, double, double, double, const dem_vec&);
    double kernel(double, double);
    void cellrange(lexer*, int, double, int&, int&, int&, int&, int&, int&);

    // output
    void print(lexer*, ghostcell*);
    void print_vtp(lexer*);
    void print_state(lexer*);
    void print_log(lexer*, ghostcell*);

    dem_core core;

    int solver;                 // 6 CFD, 5 NHFLOW
    vector<int> basemode;       // coupling mode from E 11/E 12; NHFLOW switches surface-piercing particles to unresolved
    int coupling;
    bool initialized;
    int nb;

    string infile;
    double node_spacing_factor;
    int sdf_res, quad_res;
    int gridwalls;

    // per particle fluid data
    vector<double> rhof, nuf, epsf, vsub, dvol;
    vector<dem_vec> ufl, ufl_old, Fb, Tb;
    vector<int> fluidcount;
    vector<bool> ufl_valid;

    // resolved: hydrodynamic force from the forcing, fluid mass inside the particle
    vector<dem_vec> Fibm, Tibm;
    vector<double> mfl, hvol;
    vector<dem_vec> Fstage, Tstage;     // per RK stage, combined with the RK weights at the final stage
    void combine_stages(lexer*, ghostcell*, int, double);
    bool rkwarn = false;
    double dxs = 1.0;           // horizontal grid length scale (NHFLOW: DXM mixes in the sigma spacing)
    vector<dem_vec> vprev, wprev;

    // resolved: fluid momentum and angular momentum inside the particles (Kempe & Froehlich 2012)
    vector<dem_vec> Ifl, Lfl, Ifl_old, Lfl_old;
    vector<bool> Ifl_valid;
    void internal_cfd(lexer*, fdm*, ghostcell*);
    void internal_nhflow(lexer*, fdm_nhf*, ghostcell*);
    double dt_old;

    // unresolved: reaction force on the fluid
    vector<dem_vec> Ffp;

    // hydrodynamic force on the particles over the last step (output)
    vector<dem_vec> Fhyd;

    // CFD source fields and solid fraction
    field1 *Sx;
    field2 *Sy;
    field3 *Sz;
    field4a *ALPHA;

    // NHFLOW source arrays and solid fraction
    double *SX,*SY,*SZ,*ALPHAV;

    // parameters
    double hybrid_ratio, Ca, hs_factor, travel, kernel_cells;
    int fluidacc;
    int minsub, nsub;

    // output
    double printtime;
    int printcount;
    double steptime;
    int maxiter_used;
};

#endif
