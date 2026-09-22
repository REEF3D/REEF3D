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

#include "reefmg.h"
#include "lexer.h"
#include "ghostcell.h"
#include "vec.h"
#include "matrix_diag.h"
#include "hypre_struct.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

reefmg::reefmg(lexer *p, ghostcell *pgc, int solve_input, int precon_input)
{
    //  solver_mode : 0 V-cycles only, 1 BiCGStab preconditioned by the V-cycle
    //  (default), 2 V-cycles with fallback to BiCGStab.
    //
    //  BiCGStab is the default because it costs the same as plain cycling on
    //  a well-shaped grid and degrades far more gracefully on an awkward one:
    //  on ragged local extents it needed 5 cycles where plain cycling needed 13.
    solve_mode = 1;

    presweep  = 1;
    postsweep = 1;

    //  N 12 1 restores symmetric (forward and backward) line sweeps.  The
    //  default alternates one forward pass before the coarse correction and
    //  one backward pass after, which halves the smoothing cost without
    //  costing convergence on a grid-aligned anisotropy.
    mg.set_sweepstyle(p->N12==1 ? 1 : 0);

    //  N 12 2: EXPERIMENTAL red-black line smoother with column solves batched
    //  across SIMD lanes.  Measured slower on CPU - about 1.5x per V-cycle and
    //  one extra BiCGStab iteration - because after the column factorisation
    //  the smoother is memory-bound, and red-black ordering streams the
    //  solution through memory twice per sweep.  Kept only because a GPU port
    //  needs independent columns.
    mg.set_ordering(p->N12==2 ? 1 : 0);

    //  N 14 32 stores the whole hierarchy in single precision and nothing in
    //  double: about 25% less solver memory and a faster V-cycle.  BiCGStab,
    //  the halo exchange and the convergence test stay in double, with the
    //  exact operator read from REEF3D's matrix, so the converged answer is
    //  unchanged.  Must be set before topology() calls setup().
    mg.set_precision(p->N14==32 ? 32 : 64);


    //  In fp32 mode reefmg stores no double coefficients; BiCGStab takes the
    //  exact operator from REEF3D's own matrix through fine_apply().
    Mcur=0;
    mg.set_fine_operator(this);

    p->Iarray(CVAL4, p->imax*p->jmax*(p->kmax+2));

    topology(p,pgc);
}

reefmg::~reefmg()
{
}

void reefmg::topology(lexer *p, ghostcell *pgc)
{
    int nprocs,myrank;
    MPI_Comm_size(pgc->mpi_comm,&nprocs);
    MPI_Comm_rank(pgc->mpi_comm,&myrank);

    //  The decomposition comes out of the DIVEMesh grid files, so the process
    //  topology is reconstructed from the box origins.
    int me[6];
    me[0]=p->origin_i; me[1]=p->origin_j; me[2]=p->origin_k;
    me[3]=p->knox;     me[4]=p->knoy;     me[5]=p->knoz;

    std::vector<int> all(6*nprocs,0);
    MPI_Allgather(me,6,MPI_INT,&all[0],6,MPI_INT,pgc->mpi_comm);

    std::vector<int> ox(nprocs),oy(nprocs),oz(nprocs);
    for(int r=0;r<nprocs;++r)
    {
        ox[r]=all[6*r+0]; oy[r]=all[6*r+1]; oz[r]=all[6*r+2];
    }
    std::sort(ox.begin(),ox.end()); ox.erase(std::unique(ox.begin(),ox.end()),ox.end());
    std::sort(oy.begin(),oy.end()); oy.erase(std::unique(oy.begin(),oy.end()),oy.end());
    std::sort(oz.begin(),oz.end()); oz.erase(std::unique(oz.begin(),oz.end()),oz.end());

    npx=ox.size(); npy=oy.size();
    const int npz=oz.size();

    //  Periodic boundaries wrap the matrix around the global edge; the halo
    //  exchange here has no wraparound, so the operator would silently be the
    //  wrong one.  Refuse rather than solve a different problem.
    if(p->periodic1>0 || p->periodic2>0 || p->periodic3>0)
    {
        if(p->mpirank==0)
        cout<<"REEFMG periodic boundaries are not supported yet - the halo "
            <<"exchange has no wraparound.  Use N 10 10-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2749);
    }

    //  The line relaxation solves each sigma column exactly, which requires
    //  the column to live on one rank.  For FNPF that is the natural layout.
    if(npz>1)
    {
        if(p->mpirank==0)
        cout<<"REEFMG the grid is decomposed in z ("<<npz<<" ranks).  The vertical "
            <<"line solver needs each sigma column on one rank - decompose in x and y "
            <<"only, or use N 10 10-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2750);
    }

    if(npx*npy!=nprocs)
    {
        if(p->mpirank==0)
        cout<<"REEFMG the decomposition is not a Cartesian product ("
            <<npx<<" x "<<npy<<" != "<<nprocs<<" ranks).  Use N 10 10-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2751);
    }

    cx = std::lower_bound(ox.begin(),ox.end(),p->origin_i)-ox.begin();
    cy = std::lower_bound(oy.begin(),oy.end(),p->origin_j)-oy.begin();

    //  The halo exchange sends knoy*knoz across an x face and knox*knoz across
    //  a y face, so every rank in a column must share knox and every rank in a
    //  row must share knoy.  A Cartesian decomposition gives that; check it
    //  rather than discover it as an MPI truncation error.
    {
        std::vector<int> wx(npx,-1), wy(npy,-1);
        int bad=0;

        for(int r=0;r<nprocs;++r)
        {
            const int rx=std::lower_bound(ox.begin(),ox.end(),all[6*r+0])-ox.begin();
            const int ry=std::lower_bound(oy.begin(),oy.end(),all[6*r+1])-oy.begin();

            if(wx[rx]<0) wx[rx]=all[6*r+3]; else if(wx[rx]!=all[6*r+3]) bad=1;
            if(wy[ry]<0) wy[ry]=all[6*r+4]; else if(wy[ry]!=all[6*r+4]) bad=1;
        }

        if(bad)
        {
            if(p->mpirank==0)
            cout<<"REEFMG the decomposition has inconsistent block sizes - ranks "
                <<"sharing an x position must share knox, and likewise in y.  "
                <<"Use N 10 10-19 (hypre)."<<endl;

            MPI_Abort(MPI_COMM_WORLD,-2757);
        }
    }

    //  Hand the real cell widths over: DXN[i+marge] is the width of cell i and
    //  the array carries ghost cells, so DXN+marge-1 gives i = -1 .. knox.
    //  With these the coarse operators use true agglomerate volumes and centre
    //  distances, which is what makes odd local cell numbers and stretched
    //  grids behave instead of costing iterations.
    if(!mg.setup(pgc->mpi_comm,npx,npy,cx,cy,
                 p->knox,p->knoy,p->knoz,p->gknox,p->gknoy,p->N13,
                 p->DXN+p->marge-1, p->DYN+p->marge-1))
    {
        if(p->mpirank==0)
        cout<<"REEFMG "<<mg.err()<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2752);
    }

    if(p->mpirank==0)
    {
        cout<<"REEFMG procs "<<npx<<" x "<<npy<<", local box "
            <<p->knox<<" x "<<p->knoy<<" x "<<p->knoz
            <<", multigrid levels "<<mg.levels()<<endl;

        /*if(p->knox%2!=0 || (p->knoy>1 && p->knoy%2!=0))
        cout<<"REEFMG note - odd local cell numbers give ragged coarse cells.  "
            <<"The coarse operators account for this, but even knox/knoy per rank "
            <<"still converge slightly better."<<endl;*/
    }
}

void reefmg::start(lexer *p, fdm *a, ghostcell *pgc, field &f, vec &rhsvec, int var)
{
    if(p->mpirank==0)
    cout<<"REEFMG start() not implemented - use N 10 10-19 (hypre) for this equation."<<endl;

    MPI_Abort(MPI_COMM_WORLD,-2753);
}

void reefmg::startf(lexer *p, ghostcell *pgc, field &f, vec &rhs, matrix_diag &M, int var)
{
    if(p->mpirank==0)
    cout<<"REEFMG startf() not implemented - use N 10 10-19 (hypre) for this equation."<<endl;

    MPI_Abort(MPI_COMM_WORLD,-2754);
}

void reefmg::startV(lexer *p, ghostcell *pgc, double *f, vec &rhs, matrix_diag &M, int var)
{
    //  var==44: NHFLOW potential-flow initialisation (nhflow_potential_f, I 11 1)
    if(var==44)
    start_solver44(p,pgc,f,rhs,M);

    else
    {
        if(p->mpirank==0)
        cout<<"REEFMG startV() only implemented for var==44 (NHFLOW potential)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2755);
    }
}

//  The NHFLOW potential problem differs from the Laplace and pressure
//  problems in one way that matters here: every boundary is Neumann - bed,
//  walls, and the top as a rigid lid - apart from prescribed inflow/outflow
//  fluxes.  The matrix is singular, and without a Dirichlet surface the
//  depth-uniform mode sees an unscreened horizontal Laplacian spanning the
//  whole domain.  Unlike FNPF, this problem needs the full hierarchy depth
//  and a coarsest level that genuinely resolves the smoothest mode: with a
//  shallow hierarchy BiCGStab still reaches its residual tolerance, but the
//  smooth part of the solution can be wrong, because the remaining residual
//  concentrates in exactly the modes where the error is amplified most.
//
//  So this solve runs on its own temporary hierarchy - double precision,
//  lexicographic smoothing, full depth - independent of the settings chosen
//  for the time-stepping solves, and freed afterwards.  If the decomposition
//  leaves the coarsest grid too large to resolve the depth-uniform mode, the
//  solve is handed to hypre GMRES+SMG, exactly as with N 10 1x, rather than
//  risk an inaccurate initial velocity field.  It is a one-off solve at
//  start-up, so robustness matters more than speed.
//
//  Thresholds from pure-Neumann tests at full depth and with capped depth:
//  with 64 coarsest-level sweeps accuracy failed at 512 coarsest columns;
//  with 128 sweeps it held up to 2048.  512 at 128 sweeps keeps a 4x margin.
static const long POT_COARSEST_MAX  = 512;  // columns on the coarsest level
static const int  POT_COARSE_SWEEPS = 128;

void reefmg::start_solver44(lexer *p, ghostcell *pgc, double *f, vec &rhs, matrix_diag &M)
{
    p->solveriter=0;
    starttime=pgc->timer();

    reefmg_core pot;
    pot.set_precision(64);
    pot.set_coarse_sweeps(POT_COARSE_SWEEPS);

    if(!pot.setup(pgc->mpi_comm,npx,npy,cx,cy,
                  p->knox,p->knoy,p->knoz,p->gknox,p->gknoy,0,
                  p->DXN+p->marge-1, p->DYN+p->marge-1))
    {
        if(p->mpirank==0)
        cout<<"REEFMG potential: "<<pot.err()<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2758);
    }

    const sc_level &C=pot.coarsest();

    if((long)C.gnx*C.gny > POT_COARSEST_MAX)
    {
        if(p->mpirank==0)
        cout<<"REEFMG potential (var 44): coarsest grid "<<C.gnx<<" x "<<C.gny
            <<" is too large to resolve this pure-Neumann problem reliably - "
            <<"solving it with hypre GMRES+SMG instead."<<endl;

        hypre_struct fallback(p,pgc,14,11);
        fallback.startV(p,pgc,f,rhs,M,44);
        return;
    }

    fill_matrix44(p,pot,f,rhs,M);
    pot.coarsen();

    double relres=1.0;
    p->solveriter=pot.solve_auto(p->N44,p->N46,relres,1,1,1);
    p->final_res=relres;

    fillbackvec44(p,pot,f);

    if(p->mpirank==0)
    cout<<"REEFMG potential (var 44): levels "<<pot.levels()<<"  cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres<<"  "<<setprecision(3)
        <<pgc->timer()-starttime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"REEFMG potential (var 44): WARNING - iteration limit N 46 reached.  "
        <<"If the prescribed inflow and outflow do not balance, the pure-Neumann "
        <<"system is inconsistent and cannot converge with any solver."<<endl;
}

//  Same row numbering as the assembly in nhflow_potential_f::laplace (LOOP).
//  That routine turns dry and solid cells into identity rows; they are kept
//  out of the multigrid here rather than coarsened as if they were fluid.
//  PSI is cell-centred: IJK, not FIJK as for the FNPF potential.
void reefmg::fill_matrix44(lexer *p, reefmg_core &core, double *f, vec &rhs, matrix_diag &M)
{
    sc_level &L=core.fine();

    std::fill(L.p.begin(),L.p.end(),0.0);
    std::fill(L.n.begin(),L.n.end(),0.0); std::fill(L.s.begin(),L.s.end(),0.0);
    std::fill(L.w.begin(),L.w.end(),0.0); std::fill(L.e.begin(),L.e.end(),0.0);
    std::fill(L.t.begin(),L.t.end(),0.0); std::fill(L.b.begin(),L.b.end(),0.0);
    std::fill(L.u.begin(),L.u.end(),0.0); std::fill(L.f.begin(),L.f.end(),0.0);
    std::fill(L.act.begin(),L.act.end(),0);

    count=0;
    LOOP
    {
        CVAL4[IJK]=count;
        ++count;
    }

    PLAINLOOP
    {
        const long q=L.idx(i,j,k);

        if(p->flag4[IJK]>0 && p->wet[IJ]==1 && p->DF[IJK]>0)
        {
            n=CVAL4[IJK];

            L.p[q]=M.p[n];
            L.n[q]=M.n[n];   // i+1
            L.s[q]=M.s[n];   // i-1
            L.w[q]=M.w[n];   // j+1
            L.e[q]=M.e[n];   // j-1
            L.t[q]=M.t[n];   // k+1
            L.b[q]=M.b[n];   // k-1

            L.f[q]=rhs.V[n];
            L.u[q]=f[IJK];
            L.act[q]=1;

            if(L.u[q]!=L.u[q] || L.f[q]!=L.f[q])
            p->solver_error=1;
        }
        else
        L.p[q]=1.0;      // identity row
    }
}

void reefmg::fillbackvec44(lexer *p, reefmg_core &core, double *f)
{
    sc_level &L=core.fine();

    PLAINLOOP
    PFLUIDCHECK
    f[IJK]=L.u[L.idx(i,j,k)];
}

void reefmg::startM(lexer *p, ghostcell *pgc, double *x, double *rhs, double *M, int var)
{
}

void reefmg::startF(lexer *p, ghostcell *pgc, double *f, vec &rhs, matrix_diag &M, int var)
{
    if(var==8)
    start_solver8(p,pgc,f,rhs,M,var);

    else
    {
        if(p->mpirank==0)
        cout<<"REEFMG startF() only implemented for var==8 (FNPF Laplace)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2756);
    }
}

void reefmg::start_solver8(lexer *p, ghostcell *pgc, double *f, vec &rhs,
                                     matrix_diag &M, int var)
{
    p->solveriter=0;

    starttime=pgc->timer();
    fill_matrix8(p,f,rhs,M);
    p->matrixtime+=pgc->timer()-starttime;

    //  Coarse operators are rebuilt from the fine coefficients every solve.
    //  This is a sweep over 4/3 of the fine grid with no RAP, so the moving
    //  free surface costs essentially nothing.
    starttime=pgc->timer();
    mg.coarsen();
    const double coarsentime=pgc->timer()-starttime;

    starttime=pgc->timer();
    double relres=1.0;
    p->solveriter = mg.solve_auto(p->N44,p->N46,relres,presweep,postsweep,solve_mode);
    const double solvetime=pgc->timer()-starttime;

    p->final_res=relres;

    fillbackvec8(p,f);

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"REEFMG cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres
        <<"  coarsen "<<setprecision(3)<<coarsentime
        <<" s  solve "<<setprecision(3)<<solvetime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"REEFMG WARNING - iteration limit N 46 = "<<p->N46<<" reached, res "
        <<relres<<endl;
}

void reefmg::fill_matrix8(lexer *p, double *f, vec &rhs, matrix_diag &M)
{
    sc_level &L=mg.fine();
    const bool fp32=(mg.precision()==32);

    Mcur=&M;

    //  coefficients go into exactly one precision; in fp32 mode the row map
    //  lets fine_apply() read the exact double operator from M instead
    if(fp32)
    {
        std::fill(L.pf.begin(),L.pf.end(),0.0f);
        std::fill(L.nf.begin(),L.nf.end(),0.0f); std::fill(L.sf.begin(),L.sf.end(),0.0f);
        std::fill(L.wf.begin(),L.wf.end(),0.0f); std::fill(L.ef.begin(),L.ef.end(),0.0f);
        std::fill(L.tf.begin(),L.tf.end(),0.0f); std::fill(L.bf.begin(),L.bf.end(),0.0f);

        if((long)rowmap.size()!=L.size()) rowmap.resize(L.size());
        std::fill(rowmap.begin(),rowmap.end(),-1);
    }
    else
    {
        std::fill(L.p.begin(),L.p.end(),0.0);
        std::fill(L.n.begin(),L.n.end(),0.0); std::fill(L.s.begin(),L.s.end(),0.0);
        std::fill(L.w.begin(),L.w.end(),0.0); std::fill(L.e.begin(),L.e.end(),0.0);
        std::fill(L.t.begin(),L.t.end(),0.0); std::fill(L.b.begin(),L.b.end(),0.0);
    }
    std::fill(L.u.begin(),L.u.end(),0.0); std::fill(L.f.begin(),L.f.end(),0.0);
    std::fill(L.act.begin(),L.act.end(),0);

    //  same cell numbering as the Laplace assembly in fnpf_laplace_*
    count=0;
    LOOP
    {
        CVAL4[IJK]=count;
        ++count;
    }

    PLAINLOOP
    {
        const long q=L.idx(i,j,k);

        FPWDCHECK
        {
            n=CVAL4[IJK];

            if(fp32)
            {
                L.pf[q]=(float)M.p[n];
                L.nf[q]=(float)M.n[n];   // i+1
                L.sf[q]=(float)M.s[n];   // i-1
                L.wf[q]=(float)M.w[n];   // j+1
                L.ef[q]=(float)M.e[n];   // j-1
                L.tf[q]=(float)M.t[n];   // k+1
                L.bf[q]=(float)M.b[n];   // k-1
                rowmap[q]=n;
            }
            else
            {
                L.p[q]=M.p[n];
                L.n[q]=M.n[n];   // i+1
                L.s[q]=M.s[n];   // i-1
                L.w[q]=M.w[n];   // j+1
                L.e[q]=M.e[n];   // j-1
                L.t[q]=M.t[n];   // k+1
                L.b[q]=M.b[n];   // k-1
            }

            L.f[q]=rhs.V[n];
            L.u[q]=f[FIJK];
            L.act[q]=1;

            if(L.u[q]!=L.u[q] || L.f[q]!=L.f[q])
            p->solver_error=1;
        }

        FSWDCHECK
        {
            //  identity row, passes harmlessly through the line solve
            if(fp32) L.pf[q]=1.0f;
            else     L.p [q]=1.0;
        }
    }

    //  First row of every column whose cells are all active and numbered
    //  consecutively - every fully wet column, given CVAL4's LOOP order -
    //  so that fine_apply() can run it as a plain stencil.
    if(fp32)
    {
        const int nzl=L.nz;
        colrow0.assign((long)L.nx*L.ny,-1);

        for(int ii=0;ii<L.nx;++ii)
        for(int jj=0;jj<L.ny;++jj)
        {
            const long col=L.idx(ii,jj,0);
            const int r0=rowmap[col];
            if(r0<0) continue;

            int ok=1;
            for(int kk=1;kk<nzl;++kk) if(rowmap[col+kk]!=r0+kk){ok=0; break;}
            if(ok) colrow0[(long)ii*L.ny+jj]=r0;
        }
    }
}

//  y = A x with the exact double operator, read from REEF3D's matrix_diag.
//  Used by the Krylov iteration in fp32 mode, where reefmg keeps no double
//  coefficients.  CVAL4 numbers rows in LOOP order - i, then j, then k
//  innermost - which is also reefmg's layout, so the rows are read
//  sequentially down each column.  A different row order would still be
//  correct but scatter these reads and cost several times more.
void reefmg::fine_apply(const sc_level &L,const double *x,double *y)
{
    const matrix_diag &M=*Mcur;
    const long sx=(long)(L.ny+2)*L.nz, sy=L.nz;
    const int nzl=L.nz;
    const int *rm=&rowmap[0];

    for(int ii=0;ii<L.nx;++ii)
    for(int jj=0;jj<L.ny;++jj)
    {
        const long col=L.idx(ii,jj,0);
        const int r0=colrow0[(long)ii*L.ny+jj];

        //  Fully active column with consecutive rows: contiguous slices of M,
        //  no per-cell branches, the same stencil as the internal multiply.
        //  Terms are added in the same order as the per-cell path below, so
        //  both give bit-identical results.
        if(r0>=0)
        {
            const double *P=M.p+r0, *Nn=M.n+r0, *S=M.s+r0, *W=M.w+r0, *E=M.e+r0;
            const double *T=M.t+r0, *B=M.b+r0;
            const double *xc=x+col, *xn=xc+sx, *xs=xc-sx, *xw=xc+sy, *xe=xc-sy;
            double *yc=y+col;

            for(int kk=0;kk<nzl;++kk)
            yc[kk]=P[kk]*xc[kk]+Nn[kk]*xn[kk]+S[kk]*xs[kk]+W[kk]*xw[kk]+E[kk]*xe[kk];
            for(int kk=0;kk<nzl-1;++kk) yc[kk]+=T[kk]*xc[kk+1];
            for(int kk=1;kk<nzl;  ++kk) yc[kk]+=B[kk]*xc[kk-1];
            continue;
        }

        //  general path: dry or solid cells, or non-consecutive rows
        for(int kk=0;kk<L.nz;++kk)
        {
            const long q=col+kk;
            const int r=rm[q];
            if(r<0){y[q]=x[q]; continue;}

            double v=M.p[r]*x[q]+M.n[r]*x[q+sx]+M.s[r]*x[q-sx]+M.w[r]*x[q+sy]+M.e[r]*x[q-sy];
            if(kk<L.nz-1) v+=M.t[r]*x[q+1];
            if(kk>0)      v+=M.b[r]*x[q-1];
            y[q]=v;
        }
    }
}

void reefmg::fillbackvec8(lexer *p, double *f)
{
    sc_level &L=mg.fine();

    PLAINLOOP
    {
        FPWDCHECK
        f[FIJK]=L.u[L.idx(i,j,k)];
    }
}
