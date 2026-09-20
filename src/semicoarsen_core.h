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

#ifndef SEMICOARSEN_CORE_H_
#define SEMICOARSEN_CORE_H_

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
//  benchmarked on its own; semicoarsen_fnpf.cpp is the only REEF3D-facing
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
    std::vector<double> u,f,r;
    std::vector<double> hx,hy;           // cell widths, index i+1 for i in [-1,nx]

    long idx(int i,int j,int k) const
    {
        return ((long)(i+1)*(ny+2) + (j+1))*nz + k;
    }
    long size() const {return (long)(nx+2)*(ny+2)*nz;}
};

class semicoarsen_core
{
public:

    semicoarsen_core();
    ~semicoarsen_core();

    //  world     : communicator of the run
    //  npx,npy   : process grid (the vertical must not be decomposed)
    //  cx,cy     : this rank's position in it
    //  nx,ny,nz  : local interior extent
    //  gnx,gny   : global horizontal extent
    //  Returns false with a message in err() if the layout cannot be used.
    //  dxn/dyn give the cell widths for i in [-1,nx] and j in [-1,ny], i.e.
    //  nx+2 and ny+2 entries including one ghost cell each side.  Pass NULL
    //  for a uniform grid.  The widths let the coarse operators use the real
    //  agglomerate volume and centre distance instead of assuming a factor of
    //  four, which is what makes ragged and stretched grids behave.
    bool setup(MPI_Comm world,int npx,int npy,int cx,int cy,
               int nx,int ny,int nz,int gnx,int gny,int maxlevel,
               const double *dxn=0,const double *dyn=0);

    sc_level& fine(){return lev[0];}
    int levels() const {return (int)lev.size();}
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

    void vcycle(int l,int pre,int post);
    void residual(sc_level &L);
    void halo(sc_level &L);
    double dot(const sc_level &L,const std::vector<double> &a,
                                 const std::vector<double> &b) const;

private:

    void line_gs(sc_level &L,int l,int sweeps,int dir);   // dir 0 fwd, 1 bwd, 2 symmetric
    void build_widths(const double *dxn,const double *dyn);
    void exchange_widths(sc_level &L);
    void restrict_xy(sc_level &F,sc_level &C);
    void prolong_xy(const sc_level &C,sc_level &F);
    void apply(sc_level &L,int l,const std::vector<double> &x,std::vector<double> &y);
    void precondition(const std::vector<double> &rhs,std::vector<double> &x,int pre,int post);

    std::vector<sc_level> lev;

    MPI_Comm comm;
    int myrank,nprocs;
    int nbx0,nbx1,nby0,nby1;     // neighbour ranks, MPI_PROC_NULL at the edge

    std::vector<double> sbuf,rbuf;
    std::vector<double> kr,krhat,kp,kv,ks,kt,ky,kz;   // BiCGStab work space

    int coarse_sweeps;
    int sweepstyle;          // 0: alternating fwd/bwd, 1: symmetric both ways
    int nfallback;
    char errmsg[512];
};

#endif
