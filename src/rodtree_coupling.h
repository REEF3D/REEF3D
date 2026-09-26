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
#include<string>

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

private:
    void ini_points(lexer*, ghostcell*);
    void reduce_samples(ghostcell*);
    void apply_samples();
    void finish_step(lexer*, ghostcell*, bool finalize);
    void print(lexer*);
    double kernel(double) const;
    double self_weight(double frac) const;

    rodtree rt;

    // Lagrangian points: nq per element
    std::vector<int> nq, first;                 // points per element, offset
    int npts;
    std::vector<double> buf;                    // 4 per point: u,v,w,in-fluid
    std::vector<Eigen::Vector3d> ufprev;
    double tprev;
    bool have_prev;

    double rho;
    double printtime;
    int printcount;
    std::string outdir;
};

#endif
