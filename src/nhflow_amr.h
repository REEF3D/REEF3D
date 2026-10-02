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

#ifndef NHFLOW_AMR_H_
#define NHFLOW_AMR_H_

#include"reefamr.h"
#include"nhflow_momentum_func.h"
#include"nhflow_convection.h"
#include"nhflow_timestep.h"
#include<vector>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class nhflow_signal_speed;
class nhflow_reconstruct;
class nhflow_convection;
class nhflow_diffusion;
class nhflow_pressure;
class nhflow_fsf;
class nhflow_forcing;
class nhflow_turbulence;
class vrans_nhflow;
class sixdof;
class sediment;
class ioflow;
class patchBC_interface;
class reefmg_core;

using namespace std;

//  Patch-based mesh refinement for REEF3D::NHFLOW, the NHFLOW module of REEFAMR (reefamr.h:
//  hierarchy, patch lexers, fill and exchange plans, restriction map, flux matching registry).
//
//  Level 0 is the native NHFLOW grid.  A refined patch refines x and y by 2 per level, with the
//  same sigma layers.  Each patch has its own lexer (horizontal geometry from the core, 3D flags
//  and the sigma arrays built here), its own fdm_nhf and its own instances of the NHFLOW classes
//  (momentum RK2/RK3, reconstruction, HLL/HLLC, free surface, pressure), so the kernels run
//  unchanged on every grid.
//
//  Time stepping: the level-0 momentum object hands the time step to this class (stage runner).
//  Per RK stage:
//   - the cells around the patches are filled with the stage input (coarse to fine): a patch of
//     the same level, or interpolated from the next coarser level (bicubic in x,y, layer by
//     layer; eta, U, V, W and P, then WL = eta + depth of the patch and UH = WL U ...)
//   - continuity and momentum part of the stage (phase_F, phase_M), finest patches first: the
//     flux hook of a patch records the face fluxes on its box faces (FEx/FEy, Fx/Fy of UH, VH,
//     WH, per sigma layer, and dfx/dfy), a coarse face next to a finer patch takes the mean of
//     its two fine faces -> mass and momentum conserved across the interfaces
//   - restriction of WL, eta, UH, VH, WH to the covered coarse cells (2x2 averages per layer)
//   - pressure projection on all grids together: the Poisson rows of every grid (nhflow_poisson,
//     assembled per grid), the unknowns are the leaf nodes; BiCGStab with a FAC preconditioner
//     (REEFMG V-cycle on level 0, patch-local reefmg_core), as the composite Laplace of FNPF AMR
//     (nhflow_amr_press.cpp)
//   - restriction of the corrected UH, VH, WH and P, relaxation zones and ghost cells
//   - one global time step from the finest grid (nhflow_timestep with the patch hook)
//
//  Scope of this version: static refinement boxes (A 270 levels, A 276 boxes, A 277 boxes
//  without refinement, A 275 tile width), A 510 2/3, A 511 1/2, A 514 all, A 520 0/1/2, A 512 0,
//  A 560 0, A 550 0, B 200 0, X 10 0, S 10 0, no solids (A 580 1, A 581-590), no membranes (X 330), 3D
//  grids.  Patches stay out of the relaxation zones (B 96), the in- and outflow band and dry or
//  shallow cells (they are fully wet; level 0 keeps its wetting and drying).

struct nhflow_amr_patch : public reefamr_patch
{
    fdm_nhf *d = nullptr;
    nhflow_momentum_func *pmom = nullptr;
    nhflow_signal_speed *pss = nullptr;
    nhflow_reconstruct *precon = nullptr;
    nhflow_convection *pconv = nullptr;
    nhflow_diffusion *pdiff = nullptr;
    nhflow_pressure *ppress = nullptr;
    nhflow_fsf *pfsf = nullptr;
    nhflow_forcing *pdf = nullptr;
    nhflow_turbulence *pturb = nullptr;
    vrans_nhflow *pvrans = nullptr;
    sixdof *p6dof = nullptr;
    sediment *psed = nullptr;
    ioflow *pflow = nullptr;
    patchBC_interface *pBC = nullptr;
    nhflow_stage_obj S;

    // box faces recorded by the flux hook: rec[ipol][side][r*knoz+k], r fine face index along the
    // side, k layer; side 0 low x, 1 high x, 2 low y, 3 high y; ipol 0: dfx/dfy (r only)
    vector<double> rec[5][4];

    // composite pressure
    vector<double*> kv;             // Krylov vectors, F layout
    vector<int> row;                // FIJK -> matrix row of the patch assembly, -1 none
    reefmg_core *mg = nullptr;      // patch-local multigrid of the preconditioner
    double *ptgt = nullptr;         // unknowns of the current solve (P or PCORR)
};

class nhflow_amr : public reefamr, public nhflow_stage_runner, public nhflow_flux_hook, public nhflow_timestep_hook
{
public:
    nhflow_amr(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum*, nhflow_convection*, nhflow_timestep*);
    virtual ~nhflow_amr();

    void ini(lexer*, fdm_nhf*, ghostcell*);
    void print(lexer*, fdm_nhf*, ghostcell*);
    bool active() const { return maxlev>0 && patches_total>0; }

    // the time step of all grids (nhflow_momentum_RK2/RK3::start)
    void step(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum_func*, nhflow_stage_obj&) override;

    // face fluxes of grid id (0: level 0, n: patch n-1) before the divergence
    void flux_hook(lexer*, fdm_nhf*, int, int, double*, double*) override;

    // time step: patch maxima and cell sizes
    double dt_local_max(int) override;
    void dt_cell_size(int, double&, double&) override;

    // vector space of the composite pressure solve (nhflow_amr_press.cpp, reefamr_bicgstab)
    void pr_apply(int, int);
    void pr_prec(int, int);
    double pr_dot(int, int);
    void pr_start();
    void pr_p(double, double);
    void pr_s(double);
    void pr_x(double, double);

protected:
    // REEFAMR hooks
    reefamr_patch* patch_new() override;
    void patch_objects(reefamr_patch*, ghostcell*) override;
    void patch_delete(reefamr_patch*) override;
    void tag(int, vector<unsigned char>&) override;
    void regrid_prepare(ghostcell*) override;
    void regrid_ids() override;
    void regrid_static(ghostcell*) override;
    void regrid_state(ghostcell*, vector<reefamr_patch*>&) override;
    void regrid_finish(ghostcell*, int) override;

private:
    // a patch kernel runs: ghostcell exchange off and the fdm of the ghostcell on the patch (the
    // V-type wall conditions of ghostcell read eta, dfx, UH, U ... of the fdm it holds)
    struct pscope
    {
        ghostcell *g;
        fdm_nhf *dl0;
        bool old;
        pscope(ghostcell*, fdm_nhf*, fdm_nhf*);
        ~pscope();
    };

    nhflow_amr_patch* NP(int n) { return static_cast<nhflow_amr_patch*>(P[n]); }
    static nhflow_amr_patch* NP(reefamr_patch *c) { return static_cast<nhflow_amr_patch*>(c); }

    fdm_nhf* gfd(int g) { return (g<0) ? d0 : NP(g)->d; }
    nhflow_momentum_func* gmom(int g) { return (g<0) ? mom0 : NP(g)->pmom; }
    int fidx(lexer*, int, int, int) const;      // FIJK of lexer q
    int cidx(lexer*, int, int, int) const;      // IJK of lexer q

    // patch lexer: 3D flags, boundary list, sigma arrays
    void build_lexer3D(nhflow_amr_patch&);
    void free_lexer3D(nhflow_amr_patch&);

    // interpolation from the coarser grid g: cell (ic,jc), quadrant (ox,oy)
    void pweights(lexer*, int, int, int, int, double*);
    double pq(slice&, lexer*, int, int, int, int);
    double pq3(const double*, lexer*, int, int, int, int, int, bool);
    double plin(slice&, lexer*, int, int, int, int);
    void pcol(int, int, int, int, int, const double*, int, double*);

    // stage input arrays of grid g at stage s (s<0: end of the step)
    struct stg { slice *WL; double *UH,*VH,*WH; };
    stg stage_in(int g, int s);
    stg stage_out(int g, int s);

    // fill, restriction
    void fill_stage(ghostcell*, int, int);
    void fill_bed(ghostcell*, int);
    void bc_patch(ghostcell*, nhflow_amr_patch&, int);
    void restrict_surface(int);
    void restrict_momentum(int, bool);
    void halo0(ghostcell*, int);
    void exchange_fluxes(int);
    vector<vector<double>> rval;    // [target grid id+1][rmatch index * NF]: fine face values from other ranks
    int NF = 0;                     // values per face entry: 4 variables x layers + dfx
    void prolong_patch(ghostcell*, nhflow_amr_patch&);
    template<class SEL> void fill_col(int, int, SEL);
    template<class SEL> void restrict_col(SEL);
    template<class SEL> void prolong_interior_col(nhflow_amr_patch&, SEL);
    int pord = 4;                   // prolongation: 3 biquadratic, 4 bicubic

    // composite pressure (nhflow_amr_press.cpp)
    void press_solve(lexer*, ghostcell*, int);
    double* pvec(int, int);         // grid, vector (-1: the unknowns of the solve)
    void pr_core(lexer*, ghostcell*);
    void pr_prepare(ghostcell*);
    void pr_sync(int);
    void pr_rows();
    int pr_layout = -1;
    struct pgrid
    {
        vector<int> lq, lr;         // leaf unknowns: F index, matrix row (without the fixed rows)
        vector<int> cq;             // covered unknowns: F index
        vector<int> aq, ar;         // all leaf rows
    };
    vector<pgrid> pg;               // [g+1]
    vector<double*> kv0;            // level-0 Krylov vectors
    reefmg_core *mg0 = nullptr;     // level-0 multigrid of the preconditioner
    vector<int> rowmap0;            // level 0: FIJK -> matrix row
    double *ptgt0 = nullptr;
    long pr_it_total = 0, pr_solves = 0;
    int pr_it_last = 0;
    double pr_res_last = 0.0;
    int layout_id = 0;

    // output
    void write_vtr(lexer*, nhflow_amr_patch&, int);
    void write_vtr0(lexer*, fdm_nhf*);
    void gauges(lexer*, fdm_nhf*, ghostcell*);
    double mass(lexer*, fdm_nhf*, ghostcell*);
    ofstream gaugeout, logout;

    fdm_nhf *d0;
    nhflow_momentum_func *mom0;
    nhflow_stage_obj *S0p = nullptr;
    sixdof *p6v;                    // objects shared by all patches (their NHFLOW calls do nothing)
    ioflow *pflowv;
    patchBC_interface *pBCv;
    sediment *psedv;
    nhflow_convection *pconv0;
    nhflow_timestep *pstep0;
    int cur_stage = 0;
    int gcval_eta;
    int bc5, bc6;                   // boundary types of the bed and the free surface (gcb4)
    double printtime_amr;
    int printcount_amr;
    double m0 = 0.0;
    double tm[8];
};

#endif
