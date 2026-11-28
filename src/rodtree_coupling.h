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

#ifndef RODTREE_COUPLING_H_
#define RODTREE_COUPLING_H_

// Unresolved (actuator-line type) two-way coupling of the rod-tree solver
// with REEF3D::NHFLOW and REEF3D::CFD.   ctrl.txt: Z 20 1, Z 21 dt_print
//
// Per RK stage:  sample fluid velocity at Lagrangian points (owner rank
// interpolates, one MPI_Allreduce), Morison loads on the elements, spread
// the reaction with the 4-point Peskin kernel as an acceleration.  The
// fluid-side drag is treated point-implicitly (linearised in the relative
// velocity using the kernel self-weight), so dense canopies stay stable
// at large time steps.  The structure is advanced once per time step in
// the final stage (staggered explicit coupling; the structure itself is
// implicit incl. drag and added mass).
//
// The structure is replicated on every rank (small), all ranks advance it
// with identical, all-reduced fluid data.

#include"increment.h"
#include"rodtree.h"
#include<vector>
#include<algorithm>
#include<cmath>
#include<string>
#include<functional>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class field;
class slice;

class rodtree_coupling : public increment
{
public:
    rodtree_coupling(lexer*, ghostcell*);
    virtual ~rodtree_coupling();

    // NHFLOW: applies the spread acceleration directly to U,V,W and UH,VH,WH
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, double alpha, double*, double*, double*, slice&, bool finalize);

    // CFD: adds the spread acceleration to the staggered forcing fields
    void start_cfd(lexer*, fdm*, ghostcell*, double alpha, field&, field&, field&, field&, field&, field&, bool finalize);

    // CFD mesh refinement (cfd_amr, G 7 0): level 0 runs start_cfd as above; the patches
    //  - sample the fluid velocity at the points they hold (mine: the point is a leaf point of the
    //    patch), before the momentum of level 0; start_cfd then takes these samples instead of the
    //    (restricted) velocity of level 0 under the patches
    //  - get the reaction of the elements (as evaluated by start_cfd in the same stage) on their
    //    faces, as an acceleration in Fext, Gext, Hext; width: the kernel width per direction in
    //    cells of the patch (coarser near a patch edge, as the nets)
    void amr_sample_begin();
    void amr_sample(lexer*, fdm*, const std::function<bool(const double*)> &mine);
    void amr_spread(lexer*, fdm*, const std::function<void(const double*,int*)> &width);

private:
    void ini_points(lexer*, ghostcell*);
    void reduce_samples(ghostcell*);
    void apply_samples();
    void finish_step(lexer*, ghostcell*, bool finalize);
    void print(lexer*);
    double kernel(double) const;
    double self_weight(double frac) const;
    // horizontally widened kernel for representative colonies (scale s >= 1 cells)
    double kernel_s(double r, double s) const {return kernel(r/s)/s;}
    double self_weight_s(double frac, double s) const;
    int stencil_h(double s) const {return std::max(2,(int)std::ceil(2.0*s));}

    rodtree rt;

    // Lagrangian points: nq per element
    std::vector<int> nq, first;                 // points per element, offset
    std::vector<double> wgt, shk;               // per element: representative weight, horizontal kernel scale
    int npts;
    std::vector<double> buf;                    // 4 per point: u,v,w,in-fluid
    std::vector<double> abuf;                   // mesh refinement: 5 per point, u,v,w,in-fluid,sampled by a patch
    bool amr_sampled;

    // the kernel weights of the faces in the fluid: the fraction of the full weight of the
    // stencil (D: the 1D weight functions of the three directions at the faces of one component)
    // that lies in the solids or the bed (a->solid, a->topo < 0 on both sides of the face) is given
    // to the faces in the fluid, so the reaction of an element next to a pile is not lost in it
    double fluid_scale(lexer*, fdm*, int comp, int ii, int jj, int kk, int rx, int ry, int rz,
                       const std::function<double(int,int,int,int)> &D) const;
    std::vector<Eigen::Vector3d> ufprev;
    double tprev;
    bool have_prev;

    double rho;
    double printtime;
    int printcount;
    std::string outdir;
};

#endif
