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

#include"cfd_amr.h"
#include"cfd_amr_cpm.h"
#include"reef3d_mpmd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4.h"
#include"field4a.h"
#include"turbulence.h"
#include"rans_io.h"
#include"sediment_part.h"
#include"sediment_fdm.h"
#include"CPM.h"
#include"momentum_forcing.h"
#include"reinitopo_void.h"
#include"reinitopo_RK3.h"
#include"ioflow_void.h"
#include"solver_void.h"
#include"bicgstab_ijk.h"
#include"bicgstab_ijk_2D.h"
#include"patchBC_void.h"
#include"gridfile_v2.h"
#include"dmpart_api.h"
#include"grid_helper.h"
#include"mgcslice1.h"
#include"mgcslice2.h"
#include"mgcslice4.h"
#include"reefckpt.h"
#include"ckpt_registry.h"
#include<cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<sys/stat.h>

/*--------------------------------------------------------------------
CPM on the CFD mesh refinement (validation 25: the design note; validations 29, 30).

The grid of CPM. The region of the mobile bed, Q 85 xs xe ys ye zs ze (else the bounding box of the
boxes G 15), widened to whole cells of level 0: its cells of the finest level L make a regular grid.
Its global grid file (the format of DIVEMesh M 50, gridfile_v2.h: the nodes of level L, all cells
fluid, the sides walls, periodic where level 0 is periodic and the region spans the domain) is built
in memory and decomposed over all ranks by dmpart, as REEF3D M 40 does for level 0: by the cells, or
weighted by the parcels where the bed is uneven (the bed of level 0 at the start; the parcels of the
columns at the last checkpoint at a rebalance restart, M 41 1). Every rank builds the lexer (the grid
read from memory), the ghostcell (its own communicator) and the fdm of its piece, and on them the
objects of the sediment step of a CFD run: sediment_part (CPM), the reinitialisation of the bed, an
ioflow without function, a solver (BiCGStab with the Eulerian suspension of the layer, Q 58 3, 4; else
none), and a turbulence model that holds k and its dissipation (copied).

The fluid. Every cell of the grid of CPM with its margin (global index of level L, wrapped at the
periodic sides) takes the fluid of the finest grid that holds it: a patch of level L as it is, a patch
of a coarser level or level 0 interpolated (linear in x and y on its cells and faces; vertically the
value of its cell, w linear between its faces, so that the cells of CPM at the bed take the fluid of
the first cell of that grid and not a mean with the cell in the bed; the stencils of 8 cells kept with
the plan); the items are requested once from the ranks of those grids (MPI_Alltoallv), and every step
the fluid goes over them onto the grid of CPM: u, v, w, the pressure, the level set, the density, the
viscosity, the eddy viscosity and k, its dissipation (RANS); at the start also the solids and the bed.
Where the first fluid cell of a column takes the fluid of a coarser grid, the bedload layer samples its
bed shear stress (velocity and k) at the centre of the first cell of that grid (CPM::amr_href, the least
reference height Q 66 of the columns), and so do all other columns (cfd_amr::cpm_href); with Q 86 1 the
velocity at the centre of the second cell of that grid and k of the first fluid cell of each column.

The bed. Every cell of a grid in the region gets the mean of the cells of CPM under it: the bed (topo)
of the interior cells no finer patch covers, the coarser cells under the patches by the restriction
(the mean of the children, cfd_amr::restrict_level) from level L down to level 0, the cells around the
patches filled from them; every grid then updates its direct forcing as level 0 does after a bed
change (the flags of D 22 1, gcdf_update) and the bed level of its columns (fdm::bed, read by the wall
laws of k-omega and k-epsilon with sediment). With the two-way coupling (Q 50 1, Q 11 2) the drag of the
parcels (Kc, KUx, KUy, KUz of CPM_coupling.cpp, per volume) goes onto every cell of every grid in the
region and the cell around it the same way, and the momentum of every grid takes it in its stages as
CPM::fluid_forcing does on one grid (momentum_forcing::amr_psed). Outside the region the bed is fixed
(as Q 76, Q 77 on one grid); a parcel at a side of the region that is not periodic meets a wall.

The step. At the start of the time step, where level 0 runs its sediment step: the time of level 0 onto
the lexer of CPM, the fluid onto the grid, the suspension of the layer (Q 58 3, 4: on the grid of CPM, its
top a wall: the region has to reach as high as the suspension), the sediment step of CPM, the bed (and the
drag) back. The fluid step of the refinement follows.

Adaptive refinement (G 24-27). The patches may come and go under the grid of CPM: at every regrid the
plans are made again and the bed goes onto the new patches (their bed of the grid file replaced). G 27 r
flags the cells near a mobile bed: the columns of CPM whose bed shear stress is above r times the
critical one or whose bed has moved by a cell of CPM since the start, the cells of them within two cells
of the bed.

Checkpoints (P 48, I 45): the arrays of the grid of CPM (its registry) with its piece, the scalars of its
sediment objects and run, its parcels (CPM::ckpt_rank_state); on another decomposition the cells from
the pieces that hold them, the parcels by position.

Not covered (cfd_amr::cpm_why): subcycling (G 7 1), the placement of the patches (G 40 1, 2), the hotstart
from state files (I 40), floating bodies (X 10), S 10 2. The outputs of CPM (the parcels Q 180, the log of
CPM) come from the grid of CPM; level 0 writes no bed of its own.
--------------------------------------------------------------------*/

namespace
{
// a patch kernel: ghostcell without exchange, rank-local reductions, the fdm of the patch (as cfd_amr.cpp)
struct pscope
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};

// k and its dissipation of the patches, held for the bed shear stress of the layer (Q 68)
class turbulence_amr_copy : public turbulence
{
public:
    turbulence_amr_copy(lexer *p) : kin(p), eps(p) {}
    field4 kin, eps;

    void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*) override {}
    void ktimesave(lexer*, fdm*, ghostcell*) override {}
    void etimesave(lexer*, fdm*, ghostcell*) override {}
    void isource(lexer*, fdm*) override {}
    void jsource(lexer*, fdm*) override {}
    void ksource(lexer*, fdm*) override {}
    void print_3D(lexer*, fdm*, ghostcell*, std::vector<char>&, size_t&) override {}
    void ini(lexer*, fdm*, ghostcell*) override {}
    double kinval(int i, int j, int k) override { return kin(i,j,k); }
    double epsval(int i, int j, int k) override { return eps(i,j,k); }
    void gcupdate(lexer*, fdm*, ghostcell*) override {}
    double ccipol_kinval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4_kin(kin,x,y,z); }
    double ccipol_epsval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4(eps,x,y,z); }
    double ccipol_a_kinval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4_kin(kin,x,y,z); }
    double ccipol_a_epsval(lexer *p, ghostcell*, double x, double y, double z) override { return p->ccipol4a(eps,x,y,z); }
    void kinget(int i, int j, int k, double v) override { kin(i,j,k) = v; }
    void epsget(int i, int j, int k, double v) override { eps(i,j,k) = v; }
    void name_ParaView_parallel(lexer*, ofstream&) override {}
    void name_ParaView(lexer*, std::stringstream&, int*, int&) override {}
    void offset_ParaView(lexer*, int*, int&) override {}
};

// index of a cell in the arrays of lexer q (field layout)
inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

// the cell of the arrays of grid q along direction d whose value (face: type 1, centre: 0) lies at or
// below x, and the weight of it in the linear interpolation; constant outside the arrays
void locate(const lexer *q, int d, int face, double x, int &i0, double &w0)
{
    const int m = q->margin;
    const int kn = (d==0) ? q->knox : (d==1 ? q->knoy : q->knoz);
    const double *N = (d==0) ? q->XN : (d==1 ? q->YN : q->ZN);
    const int marge = increment::marge;
    auto pos = [&](int i) { return face ? N[i+1+marge] : 0.5*(N[i+marge]+N[i+1+marge]); };

    const int lo = -m, hi = kn+m-1;     // the cells of the arrays
    if(x<=pos(lo))
    {
        i0 = lo;
        w0 = 1.0;
        return;
    }
    if(x>=pos(hi))
    {
        i0 = hi-1;
        w0 = 0.0;
        return;
    }
    int a = lo, b = hi;
    while(b-a>1)
    {
        const int c = (a+b)/2;
        if(pos(c)<=x)
        a = c;
        else
        b = c;
    }
    i0 = a;
    w0 = (pos(a+1)-x)/(pos(a+1)-pos(a));
}

// the cell of the arrays of grid q along direction d that holds x (injection), the edge cells outside the
// arrays; as a stencil of two cells (i0, i0+1) with the weight w0 of i0
void locate_cell(const lexer *q, int d, double x, int &i0, double &w0)
{
    const int m = q->margin;
    const int kn = (d==0) ? q->knox : (d==1 ? q->knoy : q->knoz);
    const double *N = (d==0) ? q->XN : (d==1 ? q->YN : q->ZN);
    const int marge = increment::marge;
    int a = -m, b = kn+m-1;
    if(x<N[a+1+marge])
    b = a;
    else if(x>=N[b+marge])
    a = b;
    else
    while(b-a>1)
    {
        const int c = (a+b)/2;
        if(N[c+marge]<=x)
        a = c;
        else
        b = c;
    }
    if(a<kn+m-1)
    {
        i0 = a;
        w0 = 1.0;
    }
    else
    {
        i0 = a-1;
        w0 = 0.0;
    }
}

// the parcels of the columns of the grid of CPM at the last checkpoint (the weights of its partition at a
// rebalance restart, M 41 1)
const char *parcel_file = "./REEF3D_CFD_AMR/REEF3D_CFD_AMR_CPM_parcels.dat";

// the work of a parcel against a cell of the grid of CPM in a step (the weights of its partition): a parcel of a
// bed at rest 0.22 cells (the flume of validation 30, CPU time with 4 times the parcels and with twice the cells),
// a moving one about three times that (validation 22), 90 % of them asleep in a bedload run
const double parcel_work = 0.25;
}

// --------------------------------------------------------------------- scope
// CPM (S 10 1, Q 10 >= 1) on the refinement: what this version covers; nullptr if it is covered
const char* cfd_amr::cpm_why(lexer *p)
{
    if(p->S10!=1 || p->Q10<1)
    return "no sediment, or CPM with the bed as topography (S 10 1, Q 10 1)";
    if(p->G15<1 && p->Q85==0)
    return "CPM: the region of the mobile bed, Q 85 xs xe ys ye zs ze or the refinement boxes G 15";
    if(p->G7!=0)
    return "CPM: one time step for all levels (G 7 0)";
    if(p->G40!=0)
    return "CPM: the patches on the ranks of the level-0 cells below (G 40 0)";
    if(p->G12>0)
    return "CPM: no refinement zone of the floating bodies (G 12)";
    if(p->I40==1)
    return "CPM: no hotstart from the state files (I 40); the checkpoints (P 48, I 45) are there";
    if(p->X10>0)
    return "CPM: no floating bodies (X 10 0)";
    return nullptr;
}

// --------------------------------------------------------------------- the grid of CPM
// the global grid file of the region (format of DIVEMesh M 50), rank 0
static std::vector<char> cpm_box_grid(lexer *p0, const int *n, const std::vector<double> *X, int px, int py, int nranks)
{
    const int xma=3, yma=3, zma=3, marge=increment::marge;
    const int NI=n[0]+2*xma, NJ=n[1]+2*yma, NK=n[2]+2*zma;
    auto c3 = [&](int i, int j, int k) { return ((size_t)(i+xma)*NJ + (size_t)(j+yma))*NK + (size_t)(k+zma); };

    // the boundary types of the sides (DIVEMesh C 11-16): walls, symmetry in y of a 2D grid; periodic
    // sides (C 21, C 22) have no boundary surfaces
    const int C11=21, C14=21, C15=21, C16=21;
    const int C12 = (n[1]==1) ? 3 : 21, C13 = (n[1]==1) ? 3 : 21;
    const int C21 = px ? 1 : 0, C22 = py ? 1 : 0, C23 = 0;

    std::vector<int> flag((size_t)NI*NJ*NK,-21), solid((size_t)NI*NJ*NK,-10), fsl((size_t)NI*NJ,-10);
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    {
        fsl[(size_t)(i+xma)*NJ + (j+yma)] = 1;
        for(int k=0; k<n[2]; ++k)
        {
            flag[c3(i,j,k)] = 1;
            solid[c3(i,j,k)] = 1;
        }
    }
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    for(int k=0; k<n[2]; ++k)
    {
        if(i==0) flag[c3(i-1,j,k)] = -C11;
        if(i==n[0]-1) flag[c3(i+1,j,k)] = -C14;
        if(j==0) flag[c3(i,j-1,k)] = -C13;
        if(j==n[1]-1) flag[c3(i,j+1,k)] = -C12;
        if(k==0) flag[c3(i,j,k-1)] = -C15;
        if(k==n[2]-1) flag[c3(i,j,k+1)] = -C16;
    }

    // the boundary surfaces (DIVEMesh surface::makesurf)
    std::vector<int> surf;
    auto add = [&](int i, int j, int k, int side, int nb)
    {
        surf.push_back(i); surf.push_back(j); surf.push_back(k); surf.push_back(side); surf.push_back(std::abs(nb));
    };
    for(int i=0; i<n[0]; ++i)
    for(int j=0; j<n[1]; ++j)
    for(int k=0; k<n[2]; ++k)
    {
        if(flag[c3(i-1,j,k)]<0 && (C21<=1 || i>0)) add(i,j,k,1,flag[c3(i-1,j,k)]);
        if(flag[c3(i+1,j,k)]<0 && (C21<=1 || i<n[0]-1)) add(i,j,k,4,flag[c3(i+1,j,k)]);
        if(flag[c3(i,j-1,k)]<0 && (C22<=1 || j>0)) add(i,j,k,3,flag[c3(i,j-1,k)]);
        if(flag[c3(i,j+1,k)]<0 && (C22<=1 || j<n[1]-1)) add(i,j,k,2,flag[c3(i,j+1,k)]);
        if(flag[c3(i,j,k-1)]<0 && (C23<=1 || k>0)) add(i,j,k,5,flag[c3(i,j,k-1)]);
        if(flag[c3(i,j,k+1)]<0 && (C23<=1 || k<n[2]-1)) add(i,j,k,6,flag[c3(i,j,k+1)]);
    }
    const int nsurf = int(surf.size()/5);

    gridv2::writer w;
    w.start(gridv2::magic_glob);

    size_t pos = w.begin("GHED");
    {
        // the decomposition: dmpart chooses (M 20 1), the cuts in x (and y in 3D)
        const int iv[] = {n[0], n[1], n[2], xma, yma, zma, marge,
                          C11, C12, C13, C14, C15, C16, C21, C22, C23,
                          p0->i_dir, p0->j_dir, p0->k_dir, 0, 0, 0,
                          1, (n[1]>1 ? 1 : 0), 0, 1, 0, 0, 0, 0, 0, nranks,
                          nsurf};
        const int nint = int(sizeof(iv)/sizeof(int));
        w.put_int(nint);
        for(int q=0; q<nint; ++q)
        w.put_int(iv[q]);

        const double dv[] = {X[0][marge+1]-X[0][marge], 0.0, 0.0, 0.0,
                             X[0][marge], X[1][marge], X[2][marge], X[0][marge+n[0]], X[1][marge+n[1]], X[2][marge+n[2]],
                             p0->global_orig_x, p0->global_orig_y, p0->alpha_grid, 0.0, 0.0, 0.0};
        const int ndbl = int(sizeof(dv)/sizeof(double));
        w.put_int(ndbl);
        for(int q=0; q<ndbl; ++q)
        w.put_double(dv[q]);
    }
    w.end(pos);

    pos = w.begin("GFLG");
    w.put_rle(flag.size(),[&](size_t q){ return flag[q]; });
    w.end(pos);

    pos = w.begin("GSOL");
    w.put_rle(solid.size(),[&](size_t q){ return solid[q]; });
    w.end(pos);

    pos = w.begin("GSLF");
    w.put_rle(fsl.size(),[&](size_t q){ return fsl[q]; });
    w.end(pos);

    pos = w.begin("NODE");
    for(int d=0; d<3; ++d)
    for(int q=0; q<n[d]+1+2*marge; ++q)
    w.put_double(X[d][q]);
    w.end(pos);

    pos = w.begin("SURF");
    w.put_table((long long)nsurf,5,[&](long long r, int c){ return surf[(size_t)r*5+c]; });
    w.end(pos);

    w.finish();
    return std::move(w.buf);
}

// the weights of the cells of the grid of CPM for its partition, (i n1 + j) n2 + k; empty: by the cells.
// A rebalance restart (M 41 1): the parcels of the columns at the last checkpoint; else the bed of
// level 0 at the start (the cells of a column below it hold Q 24 parcels each). Collective.
void cfd_amr::cpm_weights(std::vector<double> &w)
{
    cfd_amr_cpm *C = cpm;
    const int n0 = C->n[0], n1 = C->n[1], n2 = C->n[2];
    std::vector<double> col((size_t)n0*n1,0.0);     // parcels of the columns
    int ok = 0;

    if(p0->M41==1)
    {
        if(myrank==0)
        {
            std::ifstream in(parcel_file);
            int a=-1, b=-1;
            std::string line;
            std::vector<double> v;
            while(in && std::getline(in,line))
            {
                if(line.empty() || line[0]=='#')
                continue;
                std::istringstream ls(line);
                if(a<0)
                {
                    ls>>a>>b;
                    continue;
                }
                double x;
                while(ls>>x)
                v.push_back(x);
            }
            if(a==n0 && b==n1 && v.size()==(size_t)n0*n1)
            {
                col = v;
                ok = 1;
            }
        }
        MPI_Bcast(&ok,1,MPI_INT,0,reef3d_comm());
        if(ok)
        MPI_Bcast(col.data(),(int)col.size(),MPI_DOUBLE,0,reef3d_comm());
        if(myrank==0)
        cout<<"CPM on the mesh refinement: "<<(ok ? "the partition of its grid weighted by the parcels of the last checkpoint" : "no parcels of a checkpoint for its partition ("+std::string(parcel_file)+")")<<endl;
    }

    if(!ok)
    {
        // the bed of level 0: its cells in the region below the bed, per column of level 0 in the region
        int R[3], b0[3], bn[3];
        for(int d=0; d<3; ++d)
        {
            R[d] = ratio(C->L,d);
            b0[d] = C->lo[d]/R[d];
            bn[d] = C->n[d]/R[d];
        }
        std::vector<double> c0((size_t)bn[0]*bn[1],0.0);
        for(int i=0; i<p0->knox; ++i)
        for(int j=0; j<p0->knoy; ++j)
        {
            const int I = i+org[0]-b0[0], J = j+org[1]-b0[1];
            if(I<0 || I>=bn[0] || J<0 || J>=bn[1])
            continue;
            for(int k=0; k<p0->knoz; ++k)
            {
                const int K = k+org[2]-b0[2];
                if(K>=0 && K<bn[2] && a0->topo(i,j,k)<0.0)
                c0[(size_t)I*bn[1]+J] += 1.0;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,c0.data(),(int)c0.size(),MPI_DOUBLE,MPI_SUM,reef3d_comm());
        for(int i=0; i<n0; ++i)
        for(int j=0; j<n1; ++j)
        col[(size_t)i*n1+j] = c0[(size_t)(i/R[0])*bn[1]+j/R[1]]*R[2]*double(MAX(p0->Q24,1));
    }

    // the same in every column: by the cells
    double cmin = col[0], cmax = col[0];
    for(double v : col)
    {
        cmin = MIN(cmin,v);
        cmax = MAX(cmax,v);
    }
    w.clear();
    if(cmax-cmin<=1.0e-9*MAX(cmax,1.0))
    return;

    w.assign((size_t)n0*n1*n2,1.0);
    for(int i=0; i<n0; ++i)
    for(int j=0; j<n1; ++j)
    for(int k=0; k<n2; ++k)
    w[((size_t)i*n1+j)*n2+k] = 1.0 + parcel_work*col[(size_t)i*n1+j]/double(n2);
}

void cfd_amr::cpm_setup(lexer *p, fdm *a, ghostcell *pgc)
{
    cfd_amr_cpm *C = cpm = new cfd_amr_cpm;
    const int marge = increment::marge;
    C->L = maxlev;

    // the region: Q 85, else the bounding box of the boxes G 15
    double bx0[3], bx1[3];
    if(p->Q85>0)
    {
        double xs = p->Q85_xs, ys = p->Q85_ys, xe = p->Q85_xe, ye = p->Q85_ye;
        p->XYin(xs,ys);
        p->XYin(xe,ye);
        bx0[0] = MIN(xs,xe); bx1[0] = MAX(xs,xe);
        bx0[1] = MIN(ys,ye); bx1[1] = MAX(ys,ye);
        bx0[2] = p->Q85_zs; bx1[2] = p->Q85_ze;
    }
    else
    {
        bx0[0] = p->G15_xs[0]; bx1[0] = p->G15_xe[0];
        bx0[1] = p->G15_ys[0]; bx1[1] = p->G15_ye[0];
        bx0[2] = p->G15_zs[0]; bx1[2] = p->G15_ze[0];
        for(int b=1; b<p->G15; ++b)
        {
            bx0[0] = MIN(bx0[0],p->G15_xs[b]); bx1[0] = MAX(bx1[0],p->G15_xe[b]);
            bx0[1] = MIN(bx0[1],p->G15_ys[b]); bx1[1] = MAX(bx1[1],p->G15_ye[b]);
            bx0[2] = MIN(bx0[2],p->G15_zs[b]); bx1[2] = MAX(bx1[2],p->G15_ze[b]);
        }
    }

    // the level-0 cells whose centre lies in it, their cells of level L
    for(int d=0; d<3; ++d)
    {
        const int N = gn(0,d), R = ratio(C->L,d);
        int lo = N, hi = -1;
        if(rr[d]==1)
        {
            lo = 0; hi = N-1;
        }
        else
        for(int I=0; I<N; ++I)
        {
            const double xc = 0.5*(node(0,d,I)+node(0,d,I+1));
            if(xc>=bx0[d] && xc<=bx1[d])
            {
                lo = MIN(lo,I);
                hi = MAX(hi,I);
            }
        }
        if(hi<lo)
        {
            if(myrank==0)
            cout<<endl<<"!!! CPM on the mesh refinement: the region of the mobile bed (Q 85, G 15) holds no cell !!!"<<endl<<endl;
            pgc->final(true);
        }
        C->lo[d] = lo*R;
        C->n[d] = (hi-lo+1)*R;
    }

    // periodic where level 0 is periodic and the region spans the domain
    C->px = (per[0] && C->lo[0]==0 && C->n[0]==gn(C->L,0)) ? 1 : 0;
    C->py = (per[1] && C->lo[1]==0 && C->n[1]==gn(C->L,1)) ? 1 : 0;

    // the weights of its partition
    std::vector<double> wcell;
    cpm_weights(wcell);

    // the global grid file, decomposed by rank 0, the grid file of every rank sent to it
    vector<char> mine;
    {
        int ok = 1;
        if(myrank==0)
        {
            vector<double> X[3];
            for(int d=0; d<3; ++d)
            {
                X[d].resize(C->n[d]+1+2*marge);
                for(int q=0; q<(int)X[d].size(); ++q)
                X[d][q] = node(C->L,d,C->lo[d]+q-marge);
            }
            vector<char> glob = cpm_box_grid(p,C->n,X,C->px,C->py,nranks);
            vector<vector<char>> rb;
            string msg;
            if(!dmpart::build_mem(std::move(glob),nranks,rb,msg,wcell.empty() ? nullptr : &wcell))
            {
                cout<<endl<<"!!! CPM on the mesh refinement: the decomposition of its grid: "<<msg<<" !!!"<<endl<<endl;
                ok = 0;
            }
            MPI_Bcast(&ok,1,MPI_INT,0,reef3d_comm());
            if(ok)
            {
                for(int r=1; r<nranks; ++r)
                {
                    long long nb = (long long)rb[r].size();
                    MPI_Send(&nb,1,MPI_LONG_LONG,r,4731,reef3d_comm());
                    MPI_Send(rb[r].data(),int(nb),MPI_CHAR,r,4732,reef3d_comm());
                }
                mine = std::move(rb[0]);
            }
        }
        else
        {
            MPI_Bcast(&ok,1,MPI_INT,0,reef3d_comm());
            if(ok)
            {
                long long nb = 0;
                MPI_Recv(&nb,1,MPI_LONG_LONG,0,4731,reef3d_comm(),MPI_STATUS_IGNORE);
                mine.resize(size_t(nb));
                MPI_Recv(mine.data(),int(nb),MPI_CHAR,0,4732,reef3d_comm(),MPI_STATUS_IGNORE);
            }
        }
        if(!ok)
        pgc->final(true);
    }

    // the lexer, ghostcell and fdm of the grid of CPM (the steps of driver for a grid of CFD)
    lexer *q = C->p = new lexer(*p,0);
    ghostcell *g = C->pgc = new ghostcell(q);
    q->grid_inmem = true;
    q->grid_mem = std::move(mine);
    q->read_grid();
    q->vellast();
    q->lexer_ini();
    q->amr_inherit(*p);
    q->solidread = p->solidread;
    q->toporead = p->toporead;
    q->topoforcing = p->topoforcing;
    q->G501 = 1;

    g->gc_ini(q);
    g->gcx_ini(q);
    q->gridini(g);

    q->flagini();
    g->flagfield(q);
    g->flagx(q,q->flagsf1);
    g->flagx(q,q->flagsf2);
    g->flagx(q,q->flagsf3);
    g->flagx(q,q->flagsf4);
    g->flagx(q,q->flag1);
    g->flagx(q,q->flag2);
    g->flagx(q,q->flag3);
    g->flagx(q,q->flag4);
    g->gcxupdate(q);
    q->vecsize(g);
    {
        grid_helper gh(q);
        gh.fillgcb1(q);
        gh.fillgcb2(q);
        gh.fillgcb3(q);
        gh.fillgcb4a(q);
        gh.fillgcb4_wall(q);
        gh.make_dgc(q);
        gh.fill_dgc1(q);
        gh.fill_dgc2(q);
        gh.fill_dgc3(q);
        gh.fill_dgc4(q);
    }
    {
        g->gcslflagx(q,q->flagslice4);
        mgcslice1 m1(q);
        mgcslice2 m2(q);
        mgcslice4 m4(q);
        m1.makemgc(q);
        g->gcslflagx(q,q->flagslice1);
        m1.gcb_seed(q);
        m2.makemgc(q);
        g->gcslflagx(q,q->flagslice2);
        m2.gcb_seed(q);
        m4.makemgc(q);
        m4.gcb_seed(q);
        g->gcsl_setbc1(q);
        g->gcsl_setbc2(q);
        g->gcsl_setbc4(q);
        g->gcsl_setbcio(q);
        g->dgcslini1(q);
        g->dgcslini2(q);
        g->dgcslini4(q);
    }

    C->a = new fdm(q);
    g->fdm_update(C->a);

    // the sides of the cells around the piece (initialize::start): the parcels leaving it for the rank
    // across or over a periodic side go by them (part::xchange_count)
    g->flagbase(q,C->a);

    // the pieces of all ranks
    {
        int mybox[6] = {q->origin_i,q->origin_j,q->origin_k,q->knox,q->knoy,q->knoz};
        C->pbox.assign(6*nranks,0);
        MPI_Allgather(mybox,6,MPI_INT,C->pbox.data(),6,MPI_INT,reef3d_comm());
    }

    if(myrank==0)
    cout<<"CPM on the mesh refinement: the grid of level "<<C->L<<", "<<C->n[0]<<" x "<<C->n[1]<<" x "<<C->n[2]<<" cells from ("
        <<node(C->L,0,C->lo[0])<<", "<<node(C->L,1,C->lo[1])<<", "<<node(C->L,2,C->lo[2])<<") to ("
        <<node(C->L,0,C->lo[0]+C->n[0])<<", "<<node(C->L,1,C->lo[1]+C->n[1])<<", "<<node(C->L,2,C->lo[2]+C->n[2])<<")"
        <<(C->px ? ", periodic in x" : "")<<(C->py ? ", periodic in y" : "")<<(wcell.empty() ? "" : ", decomposed by the parcels")<<endl;

    // the objects of the sediment step
    C->pturb = new turbulence_amr_copy(q);
    C->pBC = new patchBC_void(q);
    C->pflow = new ioflow_v(q,g,C->pBC);
    if(p->Q58>=3)
    {
        if(q->j_dir==0)
        C->psolv = new bicgstab_ijk_2D(q,C->a,g);
        else
        C->psolv = new bicgstab_ijk(q,C->a,g);
    }
    else
    C->psolv = new solver_void(q,C->a,g);
    if(p->G540>0)
    C->preto = new reinitopo_RK3(q);
    else
    C->preto = new reinitopo_void();

    // two-way coupling: the drag of the parcels on every grid (the sediment of one grid has it in
    // momentum_forcing::psed and CPM::coupled; S 10 1 has no source of the pressure equation)
    C->coupled = (p->Q50==1 && p->Q11==2);

    // the plans, the fluid, the solids and the bed of the grids
    cpm_plans();
    cpm_bplans();
    cpm_kfields();
    cpm_pull(true);

    g->solid_forcing_ini(q,C->a);

    C->psed = new sediment_part(q,C->a,g,C->pturb,C->pBC);
    momentum_forcing::psed = nullptr;
    CPM::coupled = nullptr;
    if(C->coupled)
    momentum_forcing::amr_psed = [this](lexer *pp, fdm*, double alpha, field &u, field &v, field &w){ cpm_forcing(pp,alpha,u,v,w); };

    C->psed->ini_cfd(q,C->a,g);
    C->psed->update_cfd(q,C->a,g,C->pflow,C->preto);

    // G 27: the bed of the columns at the start (after the first sediment step)
    C->bed0.assign((size_t)q->knox*q->knoy,0.0);

    cpm_push();

    // the parcels of the ranks
    {
        vector<double> pc;
        static_cast<sediment_part*>(C->psed)->cpm_model()->parcels_per_cell(q,pc);
        double np = 0.0;
        for(double v : pc)
        np += v;
        double mx = np, sm = np;
        MPI_Allreduce(MPI_IN_PLACE,&mx,1,MPI_DOUBLE,MPI_MAX,reef3d_comm());
        MPI_Allreduce(MPI_IN_PLACE,&sm,1,MPI_DOUBLE,MPI_SUM,reef3d_comm());
        if(myrank==0 && nranks>1 && sm>0.0)
        cout<<"CPM on the mesh refinement: parcels of the ranks max/mean "<<mx/(sm/nranks)<<endl;
    }

    long nc[2] = {C->ncopy,C->nipol};
    MPI_Allreduce(MPI_IN_PLACE,nc,2,MPI_LONG,MPI_SUM,reef3d_comm());
    if(myrank==0 && nc[1]>0)
    cout<<"CPM on the mesh refinement: "<<nc[0]<<" cells of its grid (with the margin) from patches of level "<<C->L<<", "<<nc[1]<<" interpolated from coarser grids"<<endl;
}

// --------------------------------------------------------------------- the plans
// the fluid onto the grid of CPM: the finest grid at every cell of it with its margin
void cfd_amr::cpm_plans()
{
    cfd_amr_cpm *C = cpm;
    lexer *q = C->p;
    const int m = q->margin;      // the margin of the fields (3), not the nodes (marge)
    const int L = C->L;

    vector<vector<int>> req(nranks);
    vector<vector<int>> pos(nranks);
    C->ncopy = C->nipol = 0;
    C->lev_of.assign((size_t)q->knox*q->knoy*q->knoz,(signed char)L);
    int I[3], J[3];
    for(int i=-m; i<q->knox+m; ++i)
    for(int j=(q->j_dir==1 ? -m : 0); j<(q->j_dir==1 ? q->knoy+m : 1); ++j)
    for(int k=-m; k<q->knoz+m; ++k)
    {
        I[0] = C->lo[0] + q->origin_i + i;
        I[1] = C->lo[1] + q->origin_j + j;
        I[2] = C->lo[2] + q->origin_k + k;
        if(!wrap(L,I,J))
        continue;

        int gp = -1, l;
        for(l=L; l>=1; --l)
        {
            int Il[3];
            for(int d=0; d<3; ++d)
            Il[d] = J[d]/ratio(L-l,d);
            gp = gpatch_at(l,Il);
            if(gp>=0)
            break;
        }
        int rank, lid;
        if(gp>=0)
        {
            rank = GP[gp].rank;
            lid = GP[gp].lid;
        }
        else
        {
            l = 0;
            int I0[3];
            for(int d=0; d<3; ++d)
            I0[d] = J[d]/ratio(L,d);
            rank = owner0(I0);
            lid = -1;
        }
        if(l==L)
        ++C->ncopy;
        else
        ++C->nipol;
        if(i>=0 && i<q->knox && j>=0 && j<q->knoy && k>=0 && k<q->knoz)
        C->lev_of[((size_t)i*q->knoy+j)*q->knoz+k] = (signed char)l;
        const int it[5] = {lid,l,J[0],J[1],J[2]};
        req[rank].insert(req[rank].end(),it,it+5);
        pos[rank].push_back(cix(q,i,j,k));
    }

    C->nipol_all = C->nipol;
    MPI_Allreduce(MPI_IN_PLACE,&C->nipol_all,1,MPI_LONG,MPI_SUM,reef3d_comm());

    C->scnt.assign(nranks,0);
    C->rcnt.assign(nranks,0);
    C->sdsp.assign(nranks,0);
    C->rdsp.assign(nranks,0);
    for(int r=0; r<nranks; ++r)
    C->scnt[r] = (int)pos[r].size();
    MPI_Alltoall(C->scnt.data(),1,MPI_INT,C->rcnt.data(),1,MPI_INT,reef3d_comm());
    for(int r=1; r<nranks; ++r)
    {
        C->sdsp[r] = C->sdsp[r-1]+C->scnt[r-1];
        C->rdsp[r] = C->rdsp[r-1]+C->rcnt[r-1];
    }
    const int nr = C->rdsp[nranks-1]+C->rcnt[nranks-1];

    C->mine.clear();
    vector<int> sb;
    for(int r=0; r<nranks; ++r)
    {
        C->mine.insert(C->mine.end(),pos[r].begin(),pos[r].end());
        sb.insert(sb.end(),req[r].begin(),req[r].end());
    }

    vector<int> sc5(nranks), sd5(nranks), rc5(nranks), rd5(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc5[r] = 5*C->scnt[r]; sd5[r] = 5*C->sdsp[r];
        rc5[r] = 5*C->rcnt[r]; rd5[r] = 5*C->rdsp[r];
    }
    vector<int> rb(5*(size_t)nr+1);
    MPI_Alltoallv(sb.data(),sc5.data(),sd5.data(),MPI_INT,rb.data(),rc5.data(),rd5.data(),MPI_INT,reef3d_comm());

    // served: the grid and the flat index of the cell (the finest level), or the stencils of the cell
    // centre and its faces on a coarser grid
    C->sv_g.resize(nr);
    C->sv_q.resize(nr);
    C->sv_st.assign(nr,-1);
    C->st_i.clear();
    C->st_w.clear();
    for(int t=0; t<nr; ++t)
    {
        const int g = rb[5*t], l = rb[5*t+1];
        const int Jf[3] = {rb[5*t+2],rb[5*t+3],rb[5*t+4]};
        lexer *gq = glex(g);
        int og[3];
        goff(g,og);
        C->sv_g[t] = g;

        if(l==L)
        {
            C->sv_q[t] = cix(gq,Jf[0]-og[0],Jf[1]-og[1],Jf[2]-og[2]);
            continue;
        }

        // the positions of the centre and the faces of the cell of level L
        double xc[3], xf[3];
        for(int d=0; d<3; ++d)
        {
            xc[d] = 0.5*(node(L,d,Jf[d])+node(L,d,Jf[d]+1));
            xf[d] = node(L,d,Jf[d]+1);
        }
        C->sv_q[t] = -1;
        C->sv_st[t] = (int)(C->st_i.size()/32);
        for(int kind=0; kind<4; ++kind)      // u, v, w, cell centred
        {
            int i0[3];
            double w0[3];
            for(int d=0; d<3; ++d)
            {
                if(d==1 && gq->j_dir==0)
                {
                    i0[d] = 0;
                    w0[d] = 1.0;
                    continue;
                }
                // vertically the value of the coarser cell (injection; w linear between its faces): the cells of CPM
                // at the bed take the fluid of the first cell of that grid, not a mean with the cell in the bed
                if(d==2 && kind!=2)
                locate_cell(gq,d,xc[d],i0[d],w0[d]);
                else
                locate(gq,d,kind==d ? 1 : 0,kind==d ? xf[d] : xc[d],i0[d],w0[d]);
            }
            for(int a=0; a<2; ++a)
            for(int b=0; b<2; ++b)
            for(int c=0; c<2; ++c)
            {
                const int jj = (gq->j_dir==0) ? 0 : i0[1]+b;
                const double wy = (gq->j_dir==0) ? (b==0 ? 1.0 : 0.0) : (b==0 ? w0[1] : 1.0-w0[1]);
                C->st_i.push_back(cix(gq,i0[0]+a,jj,i0[2]+c));
                C->st_w.push_back((a==0 ? w0[0] : 1.0-w0[0])*wy*(c==0 ? w0[2] : 1.0-w0[2]));
            }
        }
    }
}

// the bed and the drag onto the grids: every cell of a grid in the region, the cells of CPM under it on
// the ranks that hold them
void cfd_amr::cpm_bplans()
{
    cfd_amr_cpm *C = cpm;
    const int L = C->L;
    const int mk = C->coupled ? 1 : 0;      // the drag: the cell around every grid too

    C->tg.clear(); C->tq.clear(); C->tleaf.clear(); C->tcnt.clear();
    vector<vector<int>> req(nranks), tgt(nranks);

    auto grid = [&](int g)
    {
        lexer *gq = glex(g);
        const int l = glevel(g);
        int og[3], R[3], bl[3], bh[3];
        goff(g,og);
        for(int d=0; d<3; ++d)
        {
            R[d] = ratio(L-l,d);
            bl[d] = C->lo[d]/R[d];
            bh[d] = (C->lo[d]+C->n[d])/R[d];
        }
        int I[3], J[3];
        for(int i=-mk; i<gq->knox+mk; ++i)
        for(int j=(gq->j_dir==1 ? -mk : 0); j<(gq->j_dir==1 ? gq->knoy+mk : 1); ++j)
        for(int k=-mk; k<gq->knoz+mk; ++k)
        {
            I[0] = i+og[0]; I[1] = j+og[1]; I[2] = k+og[2];
            if(!wrap(l,I,J))
            continue;
            bool in = true;
            for(int d=0; d<3; ++d)
            if(J[d]<bl[d] || J[d]>=bh[d])
            in = false;
            if(!in)
            continue;
            const bool inner = (i>=0 && i<gq->knox && j>=0 && j<gq->knoy && k>=0 && k<gq->knoz);
            const bool leaf = inner && !covered(l,J);

            // the cells of CPM under it (global index of the grid of CPM)
            int cl[3], ch[3];
            for(int d=0; d<3; ++d)
            {
                cl[d] = J[d]*R[d]-C->lo[d];
                ch[d] = cl[d]+R[d]-1;
            }
            const int t = (int)C->tg.size();
            C->tg.push_back(g);
            C->tq.push_back(cix(gq,i,j,k));
            C->tleaf.push_back(leaf ? 1 : 0);
            C->tcnt.push_back(double(R[0])*R[1]*R[2]);

            for(int r=0; r<nranks; ++r)
            {
                const int *B = &C->pbox[6*r];
                int a[3], b[3];
                bool hit = true;
                for(int d=0; d<3; ++d)
                {
                    a[d] = MAX(cl[d],B[d]);
                    b[d] = MIN(ch[d],B[d]+B[3+d]-1);
                    if(a[d]>b[d])
                    hit = false;
                }
                if(!hit)
                continue;
                const int it[6] = {a[0]-B[0],a[1]-B[1],a[2]-B[2],b[0]-B[0],b[1]-B[1],b[2]-B[2]};
                req[r].insert(req[r].end(),it,it+6);
                tgt[r].push_back(t);
            }
        }
    };
    grid(-1);
    for(int id=0; id<npat; ++id)
    grid(id);

    C->bscnt.assign(nranks,0);
    C->brcnt.assign(nranks,0);
    C->bsdsp.assign(nranks,0);
    C->brdsp.assign(nranks,0);
    for(int r=0; r<nranks; ++r)
    C->bscnt[r] = (int)tgt[r].size();
    MPI_Alltoall(C->bscnt.data(),1,MPI_INT,C->brcnt.data(),1,MPI_INT,reef3d_comm());
    for(int r=1; r<nranks; ++r)
    {
        C->bsdsp[r] = C->bsdsp[r-1]+C->bscnt[r-1];
        C->brdsp[r] = C->brdsp[r-1]+C->brcnt[r-1];
    }
    const int nr = C->brdsp[nranks-1]+C->brcnt[nranks-1];

    C->bt.clear();
    vector<int> sb;
    for(int r=0; r<nranks; ++r)
    {
        C->bt.insert(C->bt.end(),tgt[r].begin(),tgt[r].end());
        sb.insert(sb.end(),req[r].begin(),req[r].end());
    }
    vector<int> sc6(nranks), sd6(nranks), rc6(nranks), rd6(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc6[r] = 6*C->bscnt[r]; sd6[r] = 6*C->bsdsp[r];
        rc6[r] = 6*C->brcnt[r]; rd6[r] = 6*C->brdsp[r];
    }
    C->bsv.assign(6*(size_t)nr+1,0);
    MPI_Alltoallv(sb.data(),sc6.data(),sd6.data(),MPI_INT,C->bsv.data(),rc6.data(),rd6.data(),MPI_INT,reef3d_comm());
}

// two-way coupling: the drag of the parcels on every grid of this rank (zero outside the region)
void cfd_amr::cpm_kfields()
{
    cfd_amr_cpm *C = cpm;
    C->kf.clear();
    if(!C->coupled)
    return;
    auto add = [&](lexer *gq)
    {
        auto &K = C->kf[gq];
        for(int m=0; m<4; ++m)
        K[m].assign((size_t)gq->imax*gq->jmax*gq->kmax,0.0);
    };
    add(p0);
    for(int id=0; id<npat; ++id)
    add(CP(id)->pp);
}

// --------------------------------------------------------------------- the copies
// the fluid of the grids onto the grid of CPM (ini: also the solids and the bed)
void cfd_amr::cpm_pull(bool ini)
{
    cfd_amr_cpm *C = cpm;
    fdm *a = C->a;
    turbulence_amr_copy *T = static_cast<turbulence_amr_copy*>(C->pturb);

    // the fields with their kind (0, 1, 2 the faces of u, v, w, 3 cell centred)
    vector<std::function<field&(int)>> src;
    vector<field*> dst;
    vector<int> kind;
    auto add = [&](std::function<field&(int)> s, field &d, int kd) { src.push_back(s); dst.push_back(&d); kind.push_back(kd); };
    add([this](int g) -> field& { return gfd(g)->u; }, a->u, 0);
    add([this](int g) -> field& { return gfd(g)->v; }, a->v, 1);
    add([this](int g) -> field& { return gfd(g)->w; }, a->w, 2);
    add([this](int g) -> field& { return gfd(g)->press; }, a->press, 3);
    add([this](int g) -> field& { return gfd(g)->phi; }, a->phi, 3);
    add([this](int g) -> field& { return gfd(g)->ro; }, a->ro, 3);
    add([this](int g) -> field& { return gfd(g)->visc; }, a->visc, 3);
    add([this](int g) -> field& { return gfd(g)->eddyv; }, a->eddyv, 3);
    if(turb)
    {
        add([this](int g) -> field& { return kfield(g,0); }, T->kin, 3);
        add([this](int g) -> field& { return kfield(g,1); }, T->eps, 3);
    }
    if(ini)
    {
        add([this](int g) -> field& { return gfd(g)->topo; }, a->topo, 3);
        add([this](int g) -> field& { return gfd(g)->solid; }, a->solid, 3);
    }
    const int nf = (int)src.size();

    // the served items, field by field per rank (the order of the requests)
    const int nr = (int)C->sv_g.size();
    vector<double> sb((size_t)nr*nf+1);
    for(int t=0; t<nr; ++t)
    for(int f=0; f<nf; ++f)
    {
        const double *V = src[f](C->sv_g[t]).data();
        if(C->sv_st[t]<0)
        sb[(size_t)t*nf+f] = V[C->sv_q[t]];
        else
        {
            const size_t o = (size_t)C->sv_st[t]*32 + (size_t)kind[f]*8;
            double s = 0.0;
            for(int c=0; c<8; ++c)
            s += C->st_w[o+c]*V[C->st_i[o+c]];
            sb[(size_t)t*nf+f] = s;
        }
    }

    vector<int> sc(nranks), sd(nranks), rc(nranks), rd(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc[r] = nf*C->rcnt[r]; sd[r] = nf*C->rdsp[r];
        rc[r] = nf*C->scnt[r]; rd[r] = nf*C->sdsp[r];
    }
    const size_t nm = C->mine.size();
    vector<double> rb(nm*nf+1);
    MPI_Alltoallv(sb.data(),sc.data(),sd.data(),MPI_DOUBLE,rb.data(),rc.data(),rd.data(),MPI_DOUBLE,reef3d_comm());

    for(size_t t=0; t<nm; ++t)
    for(int f=0; f<nf; ++f)
    dst[f]->data()[C->mine[t]] = rb[t*nf+f];
}

// the bed of the grid of CPM onto the grids (the mean of its cells under every cell), restricted to level
// 0, the cells around the patches filled, the direct forcing and the bed level of all grids; with the
// two-way coupling the drag of the parcels as well
void cfd_amr::cpm_push()
{
    cfd_amr_cpm *C = cpm;
    fdm *a = C->a;
    lexer *q = C->p;

    vector<field*> fs;
    fs.push_back(&a->topo);
    if(C->coupled)
    {
        CPM *M = static_cast<sediment_part*>(C->psed)->cpm_model();
        for(int m=0; m<4; ++m)
        fs.push_back(&M->coupling_field(m));
    }
    const int nf = (int)fs.size();

    // the sums over the boxes of cells of the served items
    const int nr = C->brdsp[nranks-1]+C->brcnt[nranks-1];
    vector<double> sb((size_t)nr*nf+1);
    for(int t=0; t<nr; ++t)
    {
        const int *B = &C->bsv[6*(size_t)t];
        for(int f=0; f<nf; ++f)
        {
            const double *V = fs[f]->data();
            double s = 0.0;
            for(int i=B[0]; i<=B[3]; ++i)
            for(int j=B[1]; j<=B[4]; ++j)
            for(int k=B[2]; k<=B[5]; ++k)
            s += V[cix(q,i,j,k)];
            sb[(size_t)t*nf+f] = s;
        }
    }

    vector<int> sc(nranks), sd(nranks), rc(nranks), rd(nranks);
    for(int r=0; r<nranks; ++r)
    {
        sc[r] = nf*C->brcnt[r]; sd[r] = nf*C->brdsp[r];
        rc[r] = nf*C->bscnt[r]; rd[r] = nf*C->bsdsp[r];
    }
    const size_t nm = C->bt.size();
    vector<double> rb(nm*nf+1);
    MPI_Alltoallv(sb.data(),sc.data(),sd.data(),MPI_DOUBLE,rb.data(),rc.data(),rd.data(),MPI_DOUBLE,reef3d_comm());

    // the targets: the sum of their items, the mean
    const size_t nt = C->tg.size();
    vector<double> ts(nt*nf,0.0);
    for(size_t u=0; u<nm; ++u)
    for(int f=0; f<nf; ++f)
    ts[(size_t)C->bt[u]*nf+f] += rb[u*nf+f];

    if(C->coupled)
    for(auto &e : C->kf)
    for(int m=0; m<4; ++m)
    std::fill(e.second[m].begin(),e.second[m].end(),0.0);

    for(size_t t=0; t<nt; ++t)
    {
        const int g = C->tg[t];
        if(C->tleaf[t])
        gfd(g)->topo.data()[C->tq[t]] = ts[t*nf]/C->tcnt[t];
        if(C->coupled)
        {
            auto &K = C->kf[glex(g)];
            for(int m=0; m<4; ++m)
            K[m][C->tq[t]] = ts[t*nf+1+m]/C->tcnt[t];
        }
    }

    // restriction: level L down to level 0
    vector<std::function<field&(int)>> cells = {[this](int g) -> field& { return gfd(g)->topo; }};
    for(int l=C->L; l>=1; --l)
    restrict_level(l,cells,{});

    // level 0: its halo, the direct forcing
    pgc0->start4a(p0,a0->topo,150);
    if(p0->D22==1)
    pgc0->solid_forcing_flag_update(p0,a0);
    pgc0->gcdf_update(p0,a0);

    // the cells around the patches, level by level; the direct forcing of the patches
    vector<fspec> fsp(1);
    fsp[0].type = CELL;
    fsp[0].mode = 0;
    fsp[0].f = [this](int g) -> field& { return gfd(g)->topo; };
    for(int l=1; l<=maxlev; ++l)
    fill(l,fsp);

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope ps(pgc0,c->a,a0);
        if(p0->D22==1)
        pgc0->solid_forcing_flag_update(c->pp,c->a);
        pgc0->gcdf_update(c->pp,c->a);
    }

    // the bed level of the columns (the wall law of k-omega and k-epsilon samples the velocity above it
    // with sediment, S 10 > 0): from the bed of every grid
    cpm_bed(p0,a0);
    pgc0->gcsl_start4(p0,a0->bed,50);
    for(int id=0; id<npat; ++id)
    cpm_bed(CP(id)->pp,CP(id)->a);
}

// the bed level of the columns of a grid from its bed (topo): where it crosses zero (as CPM::bedzh_update),
// the bottom of the domain where the column has no bed; a column of a patch whose bed lies below it keeps
// its level
void cfd_amr::cpm_bed(lexer *pp, fdm *fa)
{
    const int m = increment::marge;
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    {
        bool set = false;
        double h = pp->ZN[0+m];
        if(fa->topo(i,j,0)>=0.0 && pp->origin_k==0)
        set = true;
        // a grid above the bottom of the domain: the bed may lie on its lower side, between the cell below it
        // (filled) and its first cell
        for(int k=(pp->origin_k>0 ? 0 : 1); k<pp->knoz; ++k)
        if(fa->topo(i,j,k-1)<0.0 && fa->topo(i,j,k)>=0.0)
        {
            h = -(fa->topo(i,j,k-1)*pp->DZP[k-1+m])/(fa->topo(i,j,k)-fa->topo(i,j,k-1)) + pp->ZP[k+m]-pp->DZP[k-1+m];
            set = true;
        }
        if(set)
        fa->bed(i,j) = h;
    }
}

// two-way coupling: the drag of the parcels in a stage of the momentum of grid pp (CPM::fluid_forcing)
void cfd_amr::cpm_forcing(lexer *pp, double alpha, field &u, field &v, field &w)
{
    auto e = cpm->kf.find(pp);
    if(e==cpm->kf.end())
    return;
    const vector<double> &Kc = e->second[0], &Kx = e->second[1], &Ky = e->second[2], &Kz = e->second[3];
    const double adt = alpha*pp->dt;
    const int di = pp->jmax*pp->kmax, dj = pp->kmax;

    for(int i=0; i<pp->knox-pp->ulast; ++i)
    for(int j=0; j<pp->knoy; ++j)
    for(int k=0; k<pp->knoz; ++k)
    {
        const int c = cix(pp,i,j,k);
        if(pp->flag1[c]>0)
        u.data()[c] = (u.data()[c] + adt*0.5*(Kx[c]+Kx[c+di]))/(1.0 + adt*0.5*(Kc[c]+Kc[c+di]));
    }

    if(pp->j_dir==1)
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy-pp->vlast; ++j)
    for(int k=0; k<pp->knoz; ++k)
    {
        const int c = cix(pp,i,j,k);
        if(pp->flag2[c]>0)
        v.data()[c] = (v.data()[c] + adt*0.5*(Ky[c]+Ky[c+dj]))/(1.0 + adt*0.5*(Kc[c]+Kc[c+dj]));
    }

    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    for(int k=0; k<pp->knoz-pp->wlast; ++k)
    {
        const int c = cix(pp,i,j,k);
        if(pp->flag3[c]>0)
        w.data()[c] = (w.data()[c] + adt*0.5*(Kz[c]+Kz[c+1]))/(1.0 + adt*0.5*(Kc[c]+Kc[c+1]));
    }
}

// --------------------------------------------------------------------- the step
void cfd_amr::cpm_step(lexer *p, fdm *a, ghostcell *pgc)
{
    if(cpm==nullptr)
    return;

    cfd_amr_cpm *C = cpm;
    lexer *q = C->p;

    // the time of level 0, and the wave period (the start of the sediment S 41 3); the sediment times
    // (sedtime, ...) are the grid's own here, not the run state of level 0 (amr_follow)
    q->wT = p->wT;
    q->simtime = p->simtime;
    q->dt = p->dt;
    q->dt_old = p->dt_old;
    q->count = p->count;

    cpm_pull(false);
    cpm_href();
    C->psed->start_susp(q,C->a,C->pgc,C->pflow,C->psolv);
    C->psed->start_cfd(q,C->a,C->pgc,C->pflow,C->preto,C->psolv);
    cpm_push();

    // G 27: the bed of the columns after the first sediment step
    if(!C->bed0_ok && q->sediter>0)
    {
        sediment_fdm *s = static_cast<sediment_part*>(C->psed)->sed_fdm();
        CPM *M = static_cast<sediment_part*>(C->psed)->cpm_model();
        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        C->bed0[(size_t)i*q->knoy+j] = M->zbl(s,i,j);
        C->bed0_ok = true;
    }
}

// the least reference height of the bed shear stress of the layer: where the first fluid cell of a column
// takes the fluid of a coarser grid (interpolated), the centre of the first cell of that grid above the bed,
// half its height (Q 66 in cells of CPM), and the same in all columns: the columns of the finer grids sample
// their flow at that height too, so that the bed shear stress does not jump at the patch edges (a fixed bed
// across a patch edge: 6 % more in the refined part with the first cell of every grid, 1 % with the common
// height, validation 30). Q 86 1: the centre of the second cell of that grid, 1.5 of its cells. Next to a patch
// the first coarse cell carries the mean of the fine cells over its height, lower than the log law at its centre:
// over a fixed bed the stress was 6 % below the uniform grids with the first cell, 1-2 % above them with the
// second, without a jump at the edges; over a mobile bed with bedforms it is 9-10 % above them (validation 34)
void cfd_amr::cpm_href()
{
    cfd_amr_cpm *C = cpm;
    lexer *q = C->p;
    CPM *M = static_cast<sediment_part*>(C->psed)->cpm_model();
    if(C->nipol_all==0)
    {
        M->amr_href.clear();
        return;
    }
    int rmax = 1;
    for(int i=0; i<q->knox; ++i)
    for(int j=0; j<q->knoy; ++j)
    {
        int k = 0;
        while(k<q->knoz-1 && C->a->topo(i,j,k)<0.0)
        ++k;
        rmax = MAX(rmax,ratio(C->L-C->lev_of[((size_t)i*q->knoy+j)*q->knoz+k],2));
    }
    MPI_Allreduce(MPI_IN_PLACE,&rmax,1,MPI_INT,MPI_MAX,reef3d_comm());
    if(rmax>1)
    M->amr_href.assign((size_t)q->knox*q->knoy,(p0->Q86==1 ? 1.5 : 0.5)*double(rmax));
    else
    M->amr_href.clear();
}

// after a regrid: the plans made again, the bed (and the drag) onto the new patches
void cfd_amr::cpm_regrid()
{
    if(cpm==nullptr)
    return;
    cpm_plans();
    cpm_bplans();
    cpm_kfields();
    cpm_push();
}

// --------------------------------------------------------------------- G 27
// the mobile columns of the grid of CPM (all ranks): the bed shear stress above G 27 times the critical
// one (of the bedload layer, Q 58 > 0, with its critical Shields number Q 60; else of the bed shear model),
// or the bed level of the column (CPM::zbl) moved by a cell of the grid since the start
void cfd_amr::cpm_mobile()
{
    cfd_amr_cpm *C = cpm;
    lexer *q = C->p;
    sediment_fdm *s = static_cast<sediment_part*>(C->psed)->sed_fdm();
    const int marge = increment::marge;
    const int NY = gn(C->L,1);
    C->mob.assign((size_t)gn(C->L,0)*NY,0);

    // the critical bed shear stress of the columns; before the sediment starts (S 41, S 45) that of the
    // critical Shields number of the layer (Q 60)
    const double tc0 = p0->Q60*(p0->S22-p0->W1)*fabs(p0->W22)*p0->S20;

    CPM *M = static_cast<sediment_part*>(C->psed)->cpm_model();
    for(int i=0; i<q->knox; ++i)
    for(int j=0; j<q->knoy; ++j)
    {
        const double b0 = C->bed0[(size_t)i*q->knoy+j];
        int k = 0;
        while(k<q->knoz-1 && q->ZN[k+1+marge]<=b0)
        ++k;
        const double dz = q->DZN[k+marge];
        const double tc = (p0->Q58==0 && s->tau_crit(i,j)>0.0) ? s->tau_crit(i,j) : tc0;
        const bool m = M->column_tau(q,s,i,j)>=p0->G27*tc || (C->bed0_ok && fabs(M->zbl(s,i,j)-b0)>=dz);
        if(m)
        {
            const int I = C->lo[0]+q->origin_i+i, J = C->lo[1]+q->origin_j+j;
            C->mob[(size_t)I*NY+J] = 1;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,C->mob.data(),(int)C->mob.size(),MPI_UNSIGNED_CHAR,MPI_MAX,reef3d_comm());
}

// G 27: cell (i,j,k) of grid g (level lc) lies within two of its cells of a mobile bed
bool cfd_amr::cpm_flag(int g, int i, int j, int k)
{
    cfd_amr_cpm *C = cpm;
    if(C==nullptr || C->mob.empty())
    return false;
    lexer *q = glex(g);
    fdm *a = gfd(g);
    const double dx = MIN3(q->DXN[i+marge],q->DZN[k+marge],(q->j_dir==1 ? q->DYN[j+marge] : 1.0e20));
    if(fabs(a->topo(i,j,k))>=2.0*dx)
    return false;

    const int l = glevel(g);
    int og[3];
    goff(g,og);
    const int R0 = ratio(C->L-l,0), R1 = ratio(C->L-l,1);
    const int NY = gn(C->L,1);
    const int I0 = (i+og[0])*R0, J0 = (j+og[1])*R1;
    for(int I=I0; I<I0+R0; ++I)
    for(int J=J0; J<J0+R1; ++J)
    if(I>=0 && I<gn(C->L,0) && J>=0 && J<NY && C->mob[(size_t)I*NY+J])
    return true;
    return false;
}

// --------------------------------------------------------------------- checkpoints
namespace
{
double code_of(const ckpt_registry::entry &e)
{
    return (e.kind=='C' ? 0.0 : (e.kind=='S' ? 2.0 : (e.kind=='F' ? 16.0 : 8.0))) + (e.isint ? 1.0 : 0.0) + (e.local ? 4.0 : 0.0);
}

long long count_of(const ckpt_registry::entry &e, const lexer *q)
{
    if(e.kind=='R')
    return e.n;
    const long long s = (long long)q->imax*q->jmax;
    return e.kind=='C' ? s*q->kmax : (e.kind=='F' ? s*q->kmaxF : s);
}

// the values per column of a field of kind k (C cells, F nodes: kmaxF = kmax+1, S one) on a grid of kmax
inline int nk_of(char k, int kmax)
{
    return k=='C' ? kmax : (k=='F' ? kmax+1 : 1);
}

// the run state of the lexer of CPM that its sediment objects change
template<class F>
void cpm_run_scalars(lexer *q, F &&f)
{
    f(q->sedtime); f(q->dtsed);
    f(q->slidecells); f(q->sedprinttime); f(q->printtime);
    double it = q->sediter;
    f(it);
    q->sediter = int(it);
}
}

void cfd_amr::cpm_ckpt_write(reefckpt &ck)
{
    cfd_amr_cpm *C = cpm;
    if(C==nullptr)
    return;
    lexer *q = C->p;

    // the grid and the run state of its lexer (the same on all ranks)
    vector<double> G = {double(C->L), double(C->lo[0]), double(C->lo[1]), double(C->lo[2]), double(C->n[0]), double(C->n[1]), double(C->n[2])};
    cpm_run_scalars(q,[&](double &v){ G.push_back(v); });
    G.push_back(C->bed0_ok ? 1.0 : 0.0);
    ck.blob("cfd_amr_cpm",G);

    // the arrays of the piece of this rank (ckpt_registry with the lexer of CPM), the bed of the start (G 27)
    vector<ckpt_registry::entry> e;
    ckpt_registry::visit(q,[&](const ckpt_registry::entry &x){ e.push_back(x); });
    vector<double> B = {double(q->origin_i), double(q->origin_j), double(q->origin_k), double(q->knox), double(q->knoy), double(q->knoz),
                        double(q->imin), double(q->imax), double(q->jmin), double(q->jmax), double(q->kmin), double(q->kmax)};
    B.push_back(double(e.size()));
    for(const auto &x : e)
    {
        const long long n = count_of(x,q);
        B.push_back(code_of(x));
        B.push_back(double(n));
        if(x.isint)
        for(long long m=0; m<n; ++m)
        B.push_back(double(((const int*)x.V)[m]));
        else
        B.insert(B.end(),(const double*)x.V,(const double*)x.V+n);
    }
    B.push_back(double(C->bed0.size()));
    B.insert(B.end(),C->bed0.begin(),C->bed0.end());
    ck.rank_blob("cfd_amr_cpm_grid",B);

    // the sediment objects: their scalars, the parcels of the rank
    vector<double> sv, sr;
    C->psed->ckpt_state(q,sv);
    ck.blob("cfd_amr_cpm_sed",sv);
    C->psed->ckpt_rank_state(q,sr);
    ck.rank_blob("cfd_amr_cpm_sed_rank",sr);

    // the parcels of the columns (the partition of a rebalance restart, M 41 1), rank 0
    {
        CPM *M = static_cast<sediment_part*>(C->psed)->cpm_model();
        vector<double> pc;
        M->parcels_per_cell(q,pc);
        vector<double> col((size_t)C->n[0]*C->n[1],0.0);
        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        for(int k=0; k<q->knoz; ++k)
        col[(size_t)(q->origin_i+i)*C->n[1]+q->origin_j+j] += pc[((size_t)i*q->knoy+j)*q->knoz+k];
        MPI_Reduce(myrank==0 ? MPI_IN_PLACE : col.data(),col.data(),(int)col.size(),MPI_DOUBLE,MPI_SUM,0,reef3d_comm());
        if(myrank==0)
        {
            mkdir("./REEF3D_CFD_AMR",0777);
            std::ofstream out(parcel_file);
            out<<"# REEF3D CFD AMR: the parcels of the columns of the grid of CPM at count "<<p0->count<<" (the partition of its grid at a rebalance restart, M 41 1)"<<endl;
            out<<C->n[0]<<" "<<C->n[1]<<endl;
            size_t m=0;
            for(double v : col)
            out<<v<<((++m)%10==0 ? "\n" : " ");
            out<<endl;
        }
    }
}

bool cfd_amr::cpm_ckpt_read(reefckpt &ck, std::string &msg)
{
    cfd_amr_cpm *C = cpm;
    if(C==nullptr)
    return true;
    lexer *q = C->p;
    ghostcell *g = C->pgc;

    const vector<double> G = ck.get_blob("cfd_amr_cpm");
    if(G.size()<7 || int(G[0])!=C->L || int(G[1])!=C->lo[0] || int(G[2])!=C->lo[1] || int(G[3])!=C->lo[2]
       || int(G[4])!=C->n[0] || int(G[5])!=C->n[1] || int(G[6])!=C->n[2])
    {
        msg = "the grid of CPM of the checkpoint differs (the region Q 85, G 15, the levels G 1), or there is none";
        return false;
    }
    {
        size_t k = 7;
        cpm_run_scalars(q,[&](double &v){ if(k<G.size()) v = G[k]; ++k; });
        if(k<G.size())
        C->bed0_ok = G[k]>0.5;
    }

    // the pieces of the checkpoint
    struct piece
    {
        int o[3], kn[3], imin, imax, jmin, jmax, kmin, kmax;
        vector<int> code;
        vector<vector<double>> data;
        vector<double> bed0;
    };
    vector<piece> pc;
    for(const auto &bl : ck.get_rank_blobs("cfd_amr_cpm_grid"))
    {
        size_t k=0;
        if(bl.size()<13)
        {
            msg = "the grid of CPM in the checkpoint is incomplete";
            return false;
        }
        piece c;
        for(int d=0; d<3; ++d) c.o[d] = int(bl[k++]);
        for(int d=0; d<3; ++d) c.kn[d] = int(bl[k++]);
        c.imin = int(bl[k++]); c.imax = int(bl[k++]); c.jmin = int(bl[k++]); c.jmax = int(bl[k++]); c.kmin = int(bl[k++]); c.kmax = int(bl[k++]);
        const int ne = int(bl[k++]);
        for(int r=0; r<ne; ++r)
        {
            if(k+2>bl.size())
            {
                msg = "the grid of CPM in the checkpoint is incomplete";
                return false;
            }
            c.code.push_back(int(bl[k++]));
            const size_t n = size_t(bl[k++]);
            if(k+n>bl.size())
            {
                msg = "the grid of CPM in the checkpoint is incomplete";
                return false;
            }
            c.data.emplace_back(bl.begin()+k,bl.begin()+k+n);
            k += n;
        }
        if(k<bl.size())
        {
            const size_t n = size_t(bl[k++]);
            if(k+n<=bl.size())
            c.bed0.assign(bl.begin()+k,bl.begin()+k+n);
        }
        pc.push_back(std::move(c));
    }

    vector<ckpt_registry::entry> e;
    ckpt_registry::visit(q,[&](const ckpt_registry::entry &x){ e.push_back(x); });

    // the same decomposition: every piece of the checkpoint is the piece of the same rank here
    int same = ((int)pc.size()==nranks) ? 1 : 0;
    if(same)
    {
        const piece &c = pc[myrank];
        if(c.o[0]!=q->origin_i || c.o[1]!=q->origin_j || c.o[2]!=q->origin_k || c.kn[0]!=q->knox || c.kn[1]!=q->knoy || c.kn[2]!=q->knoz
           || c.code.size()!=e.size())
        same = 0;
        for(size_t r=0; r<e.size() && same; ++r)
        if(c.code[r]!=int(code_of(e[r])) || (long long)c.data[r].size()!=count_of(e[r],q))
        same = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE,&same,1,MPI_INT,MPI_MIN,reef3d_comm());

    if(same)
    {
        const piece &c = pc[myrank];
        for(size_t r=0; r<e.size(); ++r)
        {
            const vector<double> &v = c.data[r];
            if(e[r].isint)
            for(size_t t=0; t<v.size(); ++t)
            ((int*)e[r].V)[t] = int(v[t]);
            else
            std::copy(v.begin(),v.end(),(double*)e[r].V);
        }
        if(c.bed0.size()==C->bed0.size())
        C->bed0 = c.bed0;
    }
    else
    {
        // another decomposition: the cells (I,J,K) of the grid from the pieces that hold them, the interior
        // of a piece first; the scalars of the objects from a piece; the solver arrays start again. The
        // arrays of a piece fit as far as both have them in the same order (as the patches, cfd_amr_ckpt.cpp)
        vector<char> fits(pc.size(),0);
        for(size_t t=0; t<pc.size(); ++t)
        {
            bool f = true;
            const size_t nc = std::min(pc[t].code.size(),e.size());
            for(size_t r=0; r<nc && f; ++r)
            if(pc[t].code[r]!=int(code_of(e[r])))
            f = false;
            fits[t] = f ? 1 : 0;
        }
        for(size_t r=0; r<e.size(); ++r)
        {
            if(e[r].kind=='R' && !e[r].local)
            {
                for(size_t t=0; t<pc.size(); ++t)
                if(fits[t] && r<pc[t].code.size() && (long long)pc[t].data[r].size()==e[r].n)
                {
                    std::copy(pc[t].data[r].begin(),pc[t].data[r].end(),(double*)e[r].V);
                    break;
                }
                continue;
            }
            if(e[r].kind=='R' || e[r].local)
            continue;
            const bool cell = (e[r].kind=='C' || e[r].kind=='F');
            const int node = (e[r].kind=='F') ? 1 : 0;   // the closing top plane of a node field
            const int nk = nk_of(e[r].kind,q->kmax);
            const int k0 = cell ? q->kmin : 0;

            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
            for(int kk=k0; kk<k0+nk; ++kk)
            {
                const int I = q->origin_i+ii, J = q->origin_j+jj, K = q->origin_k+kk;
                const piece *src = nullptr;
                int best = 2;
                for(size_t t=0; t<pc.size(); ++t)
                {
                    const piece &o = pc[t];
                    if(!fits[t] || r>=o.code.size())
                    continue;
                    const int oi = I-o.o[0], oj = J-o.o[1], ok = K-o.o[2];
                    if(oi<o.imin || oi>=o.imin+o.imax || oj<o.jmin || oj>=o.jmin+o.jmax)
                    continue;
                    if(cell && (ok<o.kmin || ok>=o.kmin+o.kmax+node))
                    continue;
                    const bool inner = (oi>=0 && oi<o.kn[0] && oj>=0 && oj<o.kn[1] && (!cell || (ok>=0 && ok<o.kn[2]+node)));
                    const int rank = inner ? 0 : 1;
                    if(rank<best)
                    {
                        best = rank;
                        src = &o;
                        if(rank==0)
                        break;
                    }
                }
                if(src==nullptr)
                continue;
                const int oi = I-src->o[0], oj = J-src->o[1], ok = K-src->o[2];
                const int snk = nk_of(e[r].kind,src->kmax);
                const size_t si = ((size_t)(oi-src->imin)*src->jmax + (oj-src->jmin))*snk + (cell ? ok-src->kmin : 0);
                const size_t di = ((size_t)(ii-q->imin)*q->jmax + (jj-q->jmin))*nk + (kk-k0);
                const double v = src->data[r][si];
                if(e[r].isint)
                ((int*)e[r].V)[di] = int(v);
                else
                ((double*)e[r].V)[di] = v;
            }
        }

        // the bed of the start of the columns
        for(int i=0; i<q->knox; ++i)
        for(int j=0; j<q->knoy; ++j)
        {
            const int I = q->origin_i+i, J = q->origin_j+j;
            for(const auto &o : pc)
            if(I>=o.o[0] && I<o.o[0]+o.kn[0] && J>=o.o[1] && J<o.o[1]+o.kn[1] && o.bed0.size()==(size_t)o.kn[0]*o.kn[1])
            {
                C->bed0[(size_t)i*q->knoy+j] = o.bed0[(size_t)(I-o.o[0])*o.kn[1]+(J-o.o[1])];
                break;
            }
        }
    }

    // the sediment objects
    if(!ck.has("cfd_amr_cpm_sed") || !C->psed->ckpt_restore(q,g,ck.get_blob("cfd_amr_cpm_sed")))
    {
        msg = "the checkpoint has no sediment state of the grid of CPM that fits this run (S 10)";
        return false;
    }
    const vector<vector<double>> sr = same ? vector<vector<double>>{ck.get_rank_blob("cfd_amr_cpm_sed_rank",myrank)}
                                           : ck.get_rank_blobs("cfd_amr_cpm_sed_rank");
    if(!C->psed->ckpt_rank_restore(q,g,sr,same!=0))
    {
        msg = "the parcels of the grid of CPM do not fit this run (S 10)";
        return false;
    }

    if(myrank==0)
    cout<<"CPM on the mesh refinement: restart from the checkpoint, "<<(same ? "its grid as it was" : "its grid mapped onto another decomposition")<<endl;

    // the bed onto the grids
    cpm_push();
    return true;
}

