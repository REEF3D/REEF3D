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

#ifndef SFLOW_AMR_H_
#define SFLOW_AMR_H_

#include"increment.h"
#include<vector>
#include<fstream>

class lexer;
class fdm2D;
class ghostcell;
class slice;
class ioflow;
class sixdof;
class patchBC_interface;
class sflow_HLL;
class sflow_signal_speed;
class sflow_reconstruct;
class sflow_diffusion;
class sflow_pressure;
class sflow_fsf;
class sflow_forcing;
class sflow_momentum_RK3;

using namespace std;

//  Patch-based mesh refinement for REEF3D::SFLOW (hydrostatic HLL, A 220 0).
//
//  Level 0 is the native SFLOW grid and stays untouched.  Each refined patch
//  is a small SFLOW domain of its own: its own lexer (control keys copied,
//  geometry of the patch), its own fdm2D and its own instances of the SFLOW
//  classes (HLL, reconstruction, signal speeds, sflow_eta, RK3), so the
//  kernels run unchanged on every grid.  This class only couples the grids:
//
//   - ghost cells of a patch are prolonged from its parent (well balanced:
//     the surface is interpolated, the depth follows the fine bed)
//   - after every RK3 stage the patch is restricted into its parent
//     (conservative averages of WL, UH, VH)
//   - the parent's faces on the patch boundary take the mean of the two fine
//     face fluxes, and the face depth of the bed-slope source is matched the
//     same way (sflow_HLL calls hll_hook between flux_bc and the divergence)
//   - one global time step, the stages of all grids in lockstep
//
//  A patch index i of the patch lexer maps to the fine cell i-1: the extra
//  cell 0 in x and y makes the kernels compute the low-side boundary face
//  (SLICELOOP1/2 start at the first cell's right face), it is refilled from
//  the parent every stage like a ghost cell.
//
//  v1: static patches (A 276 boxes, A 270 levels), rank-local and at least
//  3 parent cells away from walls, solids and partition edges.

struct sflow_amr_patch
{
    int id, lev, parent;            // parent: index into patches, -1 = level 0
    int ilo0,ihi0,jlo0,jhi0;        // region in level-0 cells (rank-local)
    int pi0,pi1,pj0,pj1;            // covered parent cells, parent lexer indices
    int nx,ny;                      // real fine cells

    lexer *pp;
    fdm2D *b;
    sflow_HLL *phll;
    sflow_signal_speed *pss;
    sflow_reconstruct *precon;
    sflow_diffusion *pdiff;
    sflow_pressure *ppress;
    sflow_fsf *pfsf;
    sflow_forcing *psfdf;
    sflow_momentum_RK3 *pmom;

    vector<int> children;

    // boundary face values recorded by hll_hook: rec[ipol][side][fine index]
    // side 0: low x, 1: high x, 2: low y, 3: high y; ipol 0 holds dfx/dfy
    vector<double> rec[5][4];
};

class sflow_amr : public increment
{
public:
    sflow_amr(lexer*, fdm2D*, ghostcell*, patchBC_interface*, sixdof*, sflow_HLL*, sflow_momentum_RK3*);
    virtual ~sflow_amr();

    void ini(lexer*, fdm2D*, ghostcell*);
    void step_begin(lexer*, fdm2D*, ghostcell*);
    void stage_begin(lexer*, fdm2D*, ghostcell*, int);
    void stage_end(lexer*, fdm2D*, ghostcell*, int);
    void step_end(lexer*, fdm2D*, ghostcell*);
    void timestep(lexer*, fdm2D*, ghostcell*);
    void print(lexer*, fdm2D*, ghostcell*);

    // called by sflow_HLL between flux_bc and the divergence (ipol 1-4)
    void hll_hook(lexer*, fdm2D*, int, int);

    int patches_total;

private:
    void make_patch(lexer*, fdm2D*, ghostcell*, int, int, int, int, int, int);
    void build_lexer(lexer*, sflow_amr_patch&, lexer*);
    void ini_bed(lexer*, sflow_amr_patch&);
    void ini_state(lexer*, ghostcell*, sflow_amr_patch&);
    void fill_ghosts(lexer*, sflow_amr_patch&, int);
    void restrict_patch(lexer*, sflow_amr_patch&, int);
    void prolong(sflow_amr_patch&, int, int, slice&, slice&, slice&, double&, double&, double&, int&);
    void write_vtr(lexer*, sflow_amr_patch&);
    void write_vtr0(lexer*, fdm2D*);
    double mass(lexer*, fdm2D*, ghostcell*);

    // parent access (level 0 or a patch)
    lexer* plex(sflow_amr_patch&);
    fdm2D* pfdm(sflow_amr_patch&);
    sflow_momentum_RK3* pmom(sflow_amr_patch&);
    int pwet(sflow_amr_patch&, int, int);
    void set_pwet(sflow_amr_patch&, int, int, int);

    vector<sflow_amr_patch> P;
    vector<int> order;              // patches sorted by level, coarse first
    int maxlev;

    lexer *p0;
    fdm2D *b0;
    sflow_HLL *phll0;
    sflow_momentum_RK3 *pmom0;
    patchBC_interface *pBC;
    sixdof *p6dof;
    ioflow *pflow_void;

    double m0, printtime_amr;
    int printcount_amr;
    ofstream logout;
    const double eps;
};

#endif
