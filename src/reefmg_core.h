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

#ifndef REEFMG_CORE_H_
#define REEFMG_CORE_H_

#include <mpi.h>
#include <vector>

//  Multigrid for the FNPF Laplace equation on a sigma grid.
//
//  The vertical is never coarsened and is solved exactly by a tridiagonal
//  sweep, so the sigma stretching drops out of the convergence rate.  Only
//  x and y are coarsened, which means the hierarchy never runs out of
//  levels no matter how few sigma layers there are.
//
//  Deliberately free of REEF3D types so that it can be compiled and
//  benchmarked on its own; reefmg.cpp is the only REEF3D-facing
//  part.
//
//  Layout on every level: k contiguous (matching REEF3D's FIJK), one halo
//  cell in x and y, none in z.  A column is therefore contiguous, which is
//  what the line solver wants.

struct sc_level
{
    int lid;               // level index, used for the message tags
    int nx,ny,nz;          // local interior extent
    int rx,ry;             // coarsening ratio from the next finer level (1 or 2)
    int gnx,gny;           // global horizontal extent (diagnostics only)

    std::vector<double> p,n,s,w,e,t,b;   // 7 diagonals, halo included
    std::vector<char>   act;             // 0 = identity row (dry, solid, padding)
    std::vector<double> tc,ti;           // factorised column matrices: modified
                                         // super-diagonal and reciprocal pivots
    std::vector<char>   colact;          // 1 if the column has any active cell
    std::vector<long>   zcol[2];         // active columns by red-black colour

    //  Coefficients live in exactly one of the two sets, depending on the
    //  storage precision chosen before setup: p..b and tc,ti in fp64 mode,
    //  pf..bf and tcf,tif in fp32 mode.  The other set stays empty.
    std::vector<float>  pf,nf,sf,wf,ef,tf,bf,tcf,tif;
    std::vector<double> u,f,r;
    std::vector<double> hx,hy;           // cell widths, index i+1 for i in [-1,nx]

    long idx(int i,int j,int k) const
    {
        return ((long)(i+1)*(ny+2) + (j+1))*nz + k;
    }
    long size() const {return (long)(nx+2)*(ny+2)*nz;}
};

//  Exact fine-level operator supplied by the host code.  In fp32 mode no
//  double-precision coefficients are stored anywhere in the hierarchy, so the
//  Krylov iteration - which must see the exact operator to guarantee its
//  tolerance - multiplies through this interface.  REEF3D implements it
//  straight from its own matrix_diag, so the matrix is never duplicated.
class sc_operator
{
public:
    virtual ~sc_operator(){}

    //  y = A x on the interior of the fine level; x already has a valid halo
    virtual void fine_apply(const sc_level &L,const double *x,double *y)=0;
};

class reefmg_core
{
public:

    reefmg_core();
    ~reefmg_core();

    //  cart      : Cartesian communicator of the run (ghostcell::cart()); the
    //              process grid, this rank's position and the neighbours are
    //              read from it.  Duplicated, not kept.  A single-rank
    //              communicator such as MPI_COMM_SELF needs no topology.
    //              The vertical must not be decomposed and nothing may be
    //              periodic.
    //  nx,ny,nz  : local interior extent
    //  gnx,gny   : global horizontal extent
    //  Returns false with a message in err() if the layout cannot be used.
    //  dxn/dyn give the cell widths for i in [-1,nx] and j in [-1,ny], i.e.
    //  nx+2 and ny+2 entries including one ghost cell each side.  Pass NULL
    //  for a uniform grid.  The widths let the coarse operators use the real
    //  agglomerate volume and centre distance instead of assuming a factor of
    //  four, which is what makes ragged and stretched grids behave.
    bool setup(MPI_Comm cart,
               int nx,int ny,int nz,int gnx,int gny,int maxlevel,
               const double *dxn=0,const double *dyn=0);

    sc_level& fine(){return lev[0];}
    const sc_level& coarsest() const {return lev.back();}
    int levels() const {return (int)lev.size();}

    //  Use only the first n levels of the hierarchy built by setup(), without
    //  rebuilding it: level n-1 becomes the coarsest and gets the coarsest-
    //  level sweeps.  For problems whose operator screens the smooth modes -
    //  a reaction term, as in the depth-averaged non-hydrostatic pressure -
    //  levels coarser than the screening length only cost time.  Clamped to
    //  [1,levels()]; setup() resets it to levels().  Takes effect at the next
    //  coarsen(), and all ranks must pass the same n.  With a truncated
    //  hierarchy the agglomerated coarse solve is bypassed.
    void set_active_levels(int n)
    {
        const int nl=(int)lev.size();
        nuse=(n<1?1:(n>nl?nl:n));
    }
    int active_levels() const {return nuse;}
    const char* err() const {return errmsg;}

    //  Build the coarse operators from the fine ones.  Call after the fine
    //  level coefficients have been filled, i.e. once per solve.
    void coarsen();

    //  BiCGStab preconditioned by one V-cycle.  Returns the iteration count;
    //  relres receives the achieved relative residual.
    int solve(double tol,int maxiter,double &relres,int pre,int post);

    //  Plain V-cycle iteration, for cases where the extra robustness of a
    //  Krylov wrapper is not needed.
    int solve_vcycle(double tol,int maxiter,double &relres,int pre,int post);

    //  mode 0: V-cycles only, 1: BiCGStab only, 2: V-cycles with automatic
    //  fallback to BiCGStab if the cycle stalls.  Mode 2 is the default: the
    //  cycle alone is cheaper, the wrapper is there for the awkward matrices.
    int solve_auto(double tol,int maxiter,double &relres,int pre,int post,int mode);

    int fallbacks() const {return nfallback;}
    void set_sweepstyle(int s){sweepstyle=s;}

    //  Storage precision of the coefficients, either/or:
    //    64: every level stored in double, as the host provides the matrix.
    //    32: every level stored in float only - no double coefficients
    //        anywhere, which roughly halves the coefficient memory.  The
    //        Krylov iteration still multiplies with the exact double
    //        operator, obtained from set_fine_operator(), so the converged
    //        answer is unchanged; only the preconditioner is approximated.
    //  Both must be called before setup(); fp32 without a fine operator is
    //  refused by setup().
    void set_precision(int bits){pcbits=(bits==32?32:64);}
    void set_fine_operator(sc_operator *op){fineop=op;}

    //  line-GS sweeps on the coarsest level (default 16).  Only matters when
    //  the coarsest grid is large - chiefly for pure-Neumann problems, where
    //  the coarsest level carries the depth-uniform mode across the domain.
    void set_coarse_sweeps(int n){coarse_sweeps=(n>0?n:16);}

    //  bytes held by the hierarchy and the Krylov work space
    long memory_bytes() const;

    //  Coarse-grid agglomeration.  The distributed hierarchy stops coarsening
    //  at a few cells per rank, which leaves a coarsest grid that grows with
    //  the number of ranks.  With agglomeration on, that coarsest level is
    //  gathered onto every rank and solved redundantly by a serial,
    //  full-depth reefmg sub-hierarchy: no scatter, no idle ranks, and the
    //  coarse matrix is gathered once per solve, the right-hand side once per
    //  V-cycle.  Required when the coarse correction must span the whole
    //  domain (pure-Neumann problems); optional otherwise.  Refused, and left
    //  off, if the gathered problem would be too large to hold on every rank.
    //  0 off, 1 on.  Must be called before setup().
    void set_agglomeration(int mode){aggmode=(mode==1?1:0);}
    void set_agglomeration_cycles(int n){agg_cycles=(n>0?n:1);}
    bool agglomerated() const {return agg!=0;}
    long agglomerated_cells() const {return agg_cells;}

    reefmg_core(const reefmg_core&)=delete;
    reefmg_core& operator=(const reefmg_core&)=delete;

    //  0: lexicographic line Gauss-Seidel, columns solved one after another.
    //  1: red-black (zebra) line Gauss-Seidel.  Columns of one colour are
    //     independent, so they are solved in batches with the tridiagonal
    //     recurrence running across SIMD lanes.
    void set_ordering(int o){ordering=(o==1?1:0);}
    int  precision() const {return pcbits;}

    void vcycle(int l,int pre,int post);
    void residual(sc_level &L);          // always double: Krylov and convergence test
    void residual_cycle(sc_level &L);    // inside the V-cycle, follows set_precision
    void halo(sc_level &L);
    void halo_vec(sc_level &L,std::vector<double> &v);
    double dot(const sc_level &L,const std::vector<double> &a,
                                 const std::vector<double> &b) const;

private:

    void line_gs(sc_level &L,int l,int sweeps,int dir);   // dir 0 fwd, 1 bwd, 2 symmetric
    template<class C> void line_gs_t(sc_level &L,int sweeps,int dir);
    template<class C> void point_gs_t(sc_level &L,int sweeps,int dir);   // nz==1
    template<class C> void line_zebra_t(sc_level &L,int sweeps,int dir);
    template<class C> void residual_t(sc_level &L);
    template<class T> void coarsen_t();
    template<class T> void factor_lines_t();
    void build_widths(const double *dxn,const double *dyn);
    void exchange_widths(sc_level &L);
    void restrict_xy(sc_level &F,sc_level &C);
    void prolong_xy(const sc_level &C,sc_level &F);
    void apply(sc_level &L,int l,std::vector<double> &x,std::vector<double> &y);
    void precondition(std::vector<double> &rhs,std::vector<double> &x,int pre,int post);
    void factor_lines();

    std::vector<sc_level> lev;

    MPI_Comm comm;
    int myrank,nprocs;
    int nbx0,nbx1,nby0,nby1;     // neighbour ranks, MPI_PROC_NULL at the edge

    std::vector<double> sbuf,rbuf;
    //  BiCGStab work space, fine level.  s shares kr - see solve() - so there
    //  are seven vectors here, not the eight the algorithm is usually written
    //  with.
    std::vector<double> kr,krhat,kp,kv,kt,ky,kz;

    int coarse_sweeps;
    int nuse;                // levels in use, see set_active_levels()
    int pcbits;              // 64 or 32: storage precision of the coefficients
    sc_operator *fineop;     // exact fine operator, required in fp32 mode
    int ordering;            // 0 lexicographic, 1 red-black batched
    std::vector<double> zr,zb,zti,ztc;   // transposed scratch for batched columns
    std::vector<int> pj,pja;             // prolongation tables for nz==1
    std::vector<double> pwy;
    int sweepstyle;          // 0: alternating fwd/bwd, 1: symmetric both ways
    int nfallback;
    char errmsg[512];

    //  coarse-grid agglomeration
    void setup_agglomeration(int npx,int npy);
    void gather_coarse_matrix();
    void coarse_solve_agg();

    reefmg_core *agg;                               // serial sub-hierarchy, on every rank
    int  aggmode, agg_cycles;
    long agg_cells;                                 // cells of the gathered problem
    std::vector<int> agg_ox,agg_oy,agg_nx,agg_ny;   // every rank's block of the coarsest level
    std::vector<int> agg_cnt,agg_disp;              // gather counts and offsets, in cells
    std::vector<double> agg_send,agg_recv;
};

#endif
