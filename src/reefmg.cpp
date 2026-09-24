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
#include "fdm.h"
#include "field.h"
#include "ghostcell.h"
#include "vec.h"
#include "matrix_diag.h"
#include "hypre_struct.h"

#include <iostream>
#include <iomanip>
#include <vector>

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

    //  N 15 1: gather the coarsest level onto every rank and solve it with a
    //  serial full-depth hierarchy.  Not needed for convergence when the free
    //  surface is Dirichlet - FNPF and the NHFLOW pressure - but it replaces
    //  the coarsest-level sweeps, each a halo exchange, by one all-gather per
    //  V-cycle, which may pay at high rank counts.  Must precede setup().
    mg.set_agglomeration(p->N15==1 ? 1 : 0);


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

    //  The process grid comes from ghostcell's Cartesian communicator, which
    //  gcx_cart_topology has already checked against the DIVEMesh neighbours.
    //  Its dimensionality follows the decomposition: 1 when only x is split,
    //  2 when z is not, so missing coordinates are 0.
    int nd=0;
    int dims[3]={1,1,1}, per[3]={0,0,0}, crd[3]={0,0,0};
    MPI_Cartdim_get(pgc->cart(),&nd);
    MPI_Cart_get(pgc->cart(),nd,dims,per,crd);

    npx = dims[0];
    npy = nd>1 ? dims[1] : 1;
    const int npz = nd>2 ? dims[2] : 1;

    cx = crd[0];
    cy = nd>1 ? crd[1] : 0;

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
            <<"line solver needs each vertical column on one rank - decompose in x and y "
            <<"only, or use N 10 10-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2750);
    }

    //  The halo exchange sends knoy*knoz across an x face and knox*knoz across
    //  a y face, so every rank in a column must share knox and every rank in a
    //  row must share knoy.  The Cartesian check compares neighbour ranks only,
    //  not block sizes, so check it rather than discover it as an MPI
    //  truncation error.
    {
        int me[4]={cx,cy,p->knox,p->knoy};
        std::vector<int> all(4*nprocs,0);
        MPI_Allgather(me,4,MPI_INT,&all[0],4,MPI_INT,pgc->mpi_comm);

        std::vector<int> wx(npx,-1), wy(npy,-1);
        int bad=0;

        for(int r=0;r<nprocs;++r)
        {
            const int rx=all[4*r+0];
            const int ry=all[4*r+1];

            if(wx[rx]<0) wx[rx]=all[4*r+2]; else if(wx[rx]!=all[4*r+2]) bad=1;
            if(wy[ry]<0) wy[ry]=all[4*r+3]; else if(wy[ry]!=all[4*r+3]) bad=1;
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
    if(!mg.setup(pgc->cart(),
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

        if(mg.agglomerated())
        cout<<"REEFMG coarse-grid agglomeration: "<<mg.agglomerated_cells()
            <<" cells gathered on every rank"<<endl;
        else if(p->N15==1 && npx*npy>1)
        cout<<"REEFMG coarse-grid agglomeration refused: the gathered problem ("
            <<mg.agglomerated_cells()<<" cells) is too large to hold on every rank"<<endl;

        /*if(p->knox%2!=0 || (p->knoy>1 && p->knoy%2!=0))
        cout<<"REEFMG note - odd local cell numbers give ragged coarse cells.  "
            <<"The coarse operators account for this, but even knox/knoy per rank "
            <<"still converge slightly better."<<endl;*/
    }
}

void reefmg::start(lexer *p, fdm *a, ghostcell *pgc, field &f, vec &rhsvec, int var)
{
    //  var==5: CFD pressure Poisson equation (pjm, pjm_corr; poisson_f and
    //  poisson_pcorr assemble a->M in LOOP order)
    if(var==5)
    start_solver5(p,pgc,f,rhsvec,a->M);

    //  var==44: CFD potential-flow initialisation (potential_f, I 11 1).
    //  var==4 reaches the pressure solver only from potential_water (I 11 2),
    //  the same kind of problem assembled the same way, so it takes the same
    //  path.  The diffusion and turbulence solves use var 1-4 as well, but
    //  they go to psolv, never to the pressure solver.
    else if(var==44 || var==4)
    start_solver44(p,a,pgc,&f,0,rhsvec,a->M);

    else
    {
        if(p->mpirank==0)
        cout<<"REEFMG start() only implemented for var==5 (CFD pressure) and var==4/44 "
            <<"(CFD potential) - use N 10 10-19 (hypre) for this equation."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2753);
    }
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
    start_solver44(p,0,pgc,0,f,rhs,M);

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
//  for the time-stepping solves, and freed afterwards.  With more than one
//  rank the coarsest level is agglomerated: gathered onto every rank and
//  solved by a serial full-depth hierarchy, which resolves the depth-uniform
//  mode however many ranks there are.  Only if the gathered problem is too
//  large to replicate, and the distributed coarsest grid too large to resolve
//  without it, is the solve handed to hypre GMRES+SMG, exactly as with
//  N 10 1x, rather than risk an inaccurate initial velocity field.
//
//  Thresholds from pure-Neumann tests without agglomeration: with 64
//  coarsest-level sweeps accuracy failed at 512 coarsest columns, with 128
//  it held up to 2048, so 512 at 128 sweeps keeps a 4x margin.  With
//  agglomeration every tested decomposition converged accurately.
static const long POT_COARSEST_MAX  = 512;  // columns on the coarsest level
static const int  POT_COARSE_SWEEPS = 128;

void reefmg::start_solver44(lexer *p, fdm *a, ghostcell *pgc, field *ff, double *f,
                            vec &rhs, matrix_diag &M)
{
    p->solveriter=0;
    starttime=pgc->timer();

    reefmg_core pot;
    pot.set_precision(64);
    pot.set_coarse_sweeps(POT_COARSE_SWEEPS);
    pot.set_agglomeration(1);

    if(!pot.setup(pgc->cart(),
                  p->knox,p->knoy,p->knoz,p->gknox,p->gknoy,0,
                  p->DXN+p->marge-1, p->DYN+p->marge-1))
    {
        if(p->mpirank==0)
        cout<<"REEFMG potential: "<<pot.err()<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2758);
    }

    const sc_level &C=pot.coarsest();

    if(!pot.agglomerated() && (long)C.gnx*C.gny > POT_COARSEST_MAX)
    {
        if(p->mpirank==0)
        cout<<"REEFMG potential (var 44): coarsest grid "<<C.gnx<<" x "<<C.gny
            <<" is too large to resolve this pure-Neumann problem reliably, and too "
            <<"large to agglomerate - solving it with hypre GMRES+SMG instead."<<endl;

        hypre_struct fallback(p,pgc,14,11);

        if(ff)
        fallback.start(p,a,pgc,*ff,rhs,44);
        else
        fallback.startV(p,pgc,f,rhs,M,44);
        return;
    }

    fill_matrix44(p,pot,ff,f,rhs,M);
    pot.coarsen();

    double relres=1.0;
    p->solveriter=pot.solve_auto(p->N44,p->N46,relres,1,1,1);
    p->final_res=relres;

    fillbackvec44(p,pot,ff,f);

    if(p->mpirank==0)
    cout<<"REEFMG potential (var 44): levels "<<pot.levels()
        <<(pot.agglomerated()? " + agglomerated coarse grid" : "")
        <<"  cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres<<"  "<<setprecision(3)
        <<pgc->timer()-starttime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"REEFMG potential (var 44): WARNING - iteration limit N 46 reached.  "
        <<"If the prescribed inflow and outflow do not balance, the pure-Neumann "
        <<"system is inconsistent and cannot converge with any solver."<<endl;
}

//  Same row numbering as the assembly in nhflow_potential_f::laplace and
//  potential_f::laplace (both LOOP).  The NHFLOW routine turns dry and solid
//  cells into identity rows; the CFD routine leaves the rows it excludes
//  (air with I 21 1, direct-forcing solids, flagsf4) empty, so there a zero
//  diagonal marks the cell as outside the problem.  Either way they are kept
//  out of the multigrid rather than coarsened as if they were fluid.
//  PSI is cell-centred: IJK, not FIJK as for the FNPF potential.
void reefmg::fill_matrix44(lexer *p, reefmg_core &core, field *ff, double *f, vec &rhs, matrix_diag &M)
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

        const bool act = ff ? (p->flag4[IJK]>0 && M.p[CVAL4[IJK]]!=0.0)                  // CFD
                            : (p->flag4[IJK]>0 && p->wet[IJ]==1 && p->DF[IJK]>0);     // NHFLOW
        if(act)
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
            L.u[q]=ff ? (*ff)(i,j,k) : f[IJK];
            L.act[q]=1;

            if(L.u[q]!=L.u[q] || L.f[q]!=L.f[q])
            p->solver_error=1;
        }
        else
        L.p[q]=1.0;      // identity row
    }
}

void reefmg::fillbackvec44(lexer *p, reefmg_core &core, field *ff, double *f)
{
    sc_level &L=core.fine();

    if(ff)
    {
        PLAINLOOP
        {
            const long q=L.idx(i,j,k);
            if(L.act[q])
            (*ff)(i,j,k)=L.u[q];
        }
        return;
    }

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
    const int nzl=L.nz;

    Mcur=&M;

    //  First call: zero everything once (this also leaves the halo of the
    //  coefficients, f and act at zero for good - nothing else writes it) and
    //  number the rows. CVAL4 depends on flag4 only, which is fixed.
    const bool first=!fill_ini;

    if(first)
    {
        if(fp32)
        {
            std::fill(L.pf.begin(),L.pf.end(),0.0f);
            std::fill(L.nf.begin(),L.nf.end(),0.0f); std::fill(L.sf.begin(),L.sf.end(),0.0f);
            std::fill(L.wf.begin(),L.wf.end(),0.0f); std::fill(L.ef.begin(),L.ef.end(),0.0f);
            std::fill(L.tf.begin(),L.tf.end(),0.0f); std::fill(L.bf.begin(),L.bf.end(),0.0f);
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

        fill_ini=true;
    }
    else
    {
        //  the solve leaves neighbour values in the halo of u; the original
        //  per-solve fill reset it to zero, so do the same for the halo only
        for(int ii=-1;ii<=L.nx;++ii)
        for(int jj=-1;jj<=L.ny;++jj)
        if(ii<0 || ii>=L.nx || jj<0 || jj>=L.ny)
        {
            double *uc=&L.u[L.idx(ii,jj,0)];
            for(int kk=0;kk<nzl;++kk) uc[kk]=0.0;
        }
    }

    //  rowmap and colrow0 depend on flag7 (fixed) and wet: rebuild them only
    //  when the wet/dry pattern has changed since the last build
    const size_t nsl=size_t(p->imax)*size_t(p->jmax);
    bool topo=first || wetsig.size()!=nsl;

    if(!topo)
    for(size_t q=0;q<nsl;++q)
    if(wetsig[q]!=p->wet[q]){topo=true; break;}

    if(topo)
    {
        wetsig.assign(p->wet,p->wet+nsl);

        if(fp32)
        {
            if((long)rowmap.size()!=L.size()) rowmap.resize(L.size());
            std::fill(rowmap.begin(),rowmap.end(),-1);
        }
    }

    //  One pass per column. Every interior entry is written each call, so the
    //  result is the same as zero-filling everything and writing the active
    //  rows, as before. Rows of a column are consecutive in M, f and CVAL4.
    const double *const Mp=M.p.data(), *const Mn=M.n.data(), *const Ms=M.s.data(), *const Mw=M.w.data();
    const double *const Me=M.e.data(), *const Mt=M.t.data(), *const Mb=M.b.data(), *const R=rhs.V.data();
    const int *const flag7=p->flag7;
    int err=0;

    ILOOP
    JLOOP
    {
        k=0;
        const int fc=FIJK;             // f / flag7 index of k=0
        const int cc=IJK;              // CVAL4 index of k=0
        const int w=p->wet[IJ];
        const long q0=L.idx(i,j,0);

        for(int kk=0;kk<nzl;++kk)
        {
            const long q=q0+kk;

            if(flag7[fc+kk]>0 && w>0)                 // FPWDCHECK
            {
                const int r=CVAL4[cc+kk];

                if(fp32)
                {
                    L.pf[q]=(float)Mp[r];
                    L.nf[q]=(float)Mn[r];   // i+1
                    L.sf[q]=(float)Ms[r];   // i-1
                    L.wf[q]=(float)Mw[r];   // j+1
                    L.ef[q]=(float)Me[r];   // j-1
                    L.tf[q]=(float)Mt[r];   // k+1
                    L.bf[q]=(float)Mb[r];   // k-1
                    if(topo) rowmap[q]=r;
                }
                else
                {
                    L.p[q]=Mp[r];
                    L.n[q]=Mn[r];   // i+1
                    L.s[q]=Ms[r];   // i-1
                    L.w[q]=Mw[r];   // j+1
                    L.e[q]=Me[r];   // j-1
                    L.t[q]=Mt[r];   // k+1
                    L.b[q]=Mb[r];   // k-1
                }

                const double fv=R[r], uv=f[fc+kk];
                L.f[q]=fv;
                L.u[q]=uv;
                L.act[q]=1;

                if(uv!=uv || fv!=fv)
                err=1;
            }
            else
            {
                //  FSWDCHECK: identity row, passes harmlessly through the
                //  line solve; any other cell: all zero (as after the fill)
                const bool ident=(flag7[fc+kk]<=0 || w==0);

                if(fp32)
                {
                    L.pf[q]=ident?1.0f:0.0f;
                    L.nf[q]=L.sf[q]=L.wf[q]=L.ef[q]=L.tf[q]=L.bf[q]=0.0f;
                }
                else
                {
                    L.p[q]=ident?1.0:0.0;
                    L.n[q]=L.s[q]=L.w[q]=L.e[q]=L.t[q]=L.b[q]=0.0;
                }

                L.f[q]=0.0;
                L.u[q]=0.0;
                L.act[q]=0;
            }
        }
    }

    if(err)
    p->solver_error=1;

    //  First row of every column whose cells are all active and numbered
    //  consecutively - every fully wet column, given CVAL4's LOOP order -
    //  so that fine_apply() can run it as a plain stencil.
    if(fp32 && topo)
    build_colrow0(L);
}

void reefmg::build_colrow0(const sc_level &L)
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
            const double *P=M.p.data()+r0, *Nn=M.n.data()+r0, *S=M.s.data()+r0, *W=M.w.data()+r0, *E=M.e.data()+r0;
            const double *T=M.t.data()+r0, *B=M.b.data()+r0;
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

//  ---------------------------------------------------------------------------
//  REEF3D::CFD pressure Poisson equation, start(...,5)
//
//  The operator is poisson_f / poisson_pcorr: a 7-point stencil in 1/rho on
//  the cell-centred Cartesian grid, rows in LOOP order (flag4>0).  Couplings
//  to flag4<0 neighbours - walls, the bed, embedded solids, the global edge -
//  have been moved to the right-hand side with the ghost values, while the
//  diagonal keeps them.  Every row next to such a cell is therefore strictly
//  diagonally dominant and the system is nonsingular, just like the FNPF
//  Laplace equation with its Dirichlet free surface: the solve runs on the
//  time-stepping hierarchy, no agglomeration needed.
//
//  Two things differ from FNPF:
//   - the active cells come from flag4, which changes with moving bodies
//     and topography updates, so rows are renumbered and rowmap/colrow0
//     rebuilt on every solve - one pass over the grid, as hypre_struct
//     does in fill_matrix4;
//   - a column can hold solid cells anywhere, e.g. under a variable bed.
//     They become identity rows, which the column solve passes through
//     and the coarsening leaves out (act = 0).
//  The density jump across the free surface lies mostly along the columns,
//  where the line solve treats it exactly.
//  ---------------------------------------------------------------------------

void reefmg::start_solver5(lexer *p, ghostcell *pgc, field &f, vec &rhs, matrix_diag &M)
{
    p->solveriter=0;

    starttime=pgc->timer();
    fill_matrix5(p,f,rhs,M);
    p->matrixtime+=pgc->timer()-starttime;

    starttime=pgc->timer();
    mg.coarsen();
    const double coarsentime=pgc->timer()-starttime;

    starttime=pgc->timer();
    double relres=1.0;
    p->solveriter = mg.solve_auto(p->N44,p->N46,relres,presweep,postsweep,solve_mode);
    const double solvetime=pgc->timer()-starttime;

    p->final_res=relres;

    fillbackvec5(p,f);

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"REEFMG cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres
        <<"  coarsen "<<setprecision(3)<<coarsentime
        <<" s  solve "<<setprecision(3)<<solvetime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"REEFMG WARNING - iteration limit N 46 = "<<p->N46<<" reached, res "
        <<relres<<endl;
}

void reefmg::fill_matrix5(lexer *p, field &f, vec &rhs, matrix_diag &M)
{
    sc_level &L=mg.fine();
    const bool fp32=(mg.precision()==32);

    Mcur=&M;

    //  First call: zero everything once, so the halo of the coefficients,
    //  f and act stays at zero for good.  Later calls reset the halo of u
    //  only; the interior is rewritten below.
    if(!fill_ini)
    {
        if(fp32)
        {
            std::fill(L.pf.begin(),L.pf.end(),0.0f);
            std::fill(L.nf.begin(),L.nf.end(),0.0f); std::fill(L.sf.begin(),L.sf.end(),0.0f);
            std::fill(L.wf.begin(),L.wf.end(),0.0f); std::fill(L.ef.begin(),L.ef.end(),0.0f);
            std::fill(L.tf.begin(),L.tf.end(),0.0f); std::fill(L.bf.begin(),L.bf.end(),0.0f);
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

        fill_ini=true;
    }
    else
    {
        for(int ii=-1;ii<=L.nx;++ii)
        for(int jj=-1;jj<=L.ny;++jj)
        if(ii<0 || ii>=L.nx || jj<0 || jj>=L.ny)
        {
            double *uc=&L.u[L.idx(ii,jj,0)];
            for(int kk=0;kk<L.nz;++kk) uc[kk]=0.0;
        }
    }

    //  same row numbering as poisson_f::start (LOOP), redone every solve
    //  because flag4 can change between time steps
    count=0;
    LOOP
    {
        CVAL4[IJK]=count;
        ++count;
    }

    if(fp32)
    {
        if((long)rowmap.size()!=L.size()) rowmap.resize(L.size());
        std::fill(rowmap.begin(),rowmap.end(),-1);
    }

    int err=0;

    PLAINLOOP
    {
        const long q=L.idx(i,j,k);

        //  a zero diagonal would mean a row poisson_f left empty; keep it
        //  out of the problem rather than hand the smoother a zero pivot
        int r=-1;
        if(p->flag4[IJK]>0)
        {
            r=CVAL4[IJK];
            if(M.p[r]==0.0) r=-1;
        }

        if(r>=0)
        {
            if(fp32)
            {
                L.pf[q]=(float)M.p[r];
                L.nf[q]=(float)M.n[r];   // i+1
                L.sf[q]=(float)M.s[r];   // i-1
                L.wf[q]=(float)M.w[r];   // j+1
                L.ef[q]=(float)M.e[r];   // j-1
                L.tf[q]=(float)M.t[r];   // k+1
                L.bf[q]=(float)M.b[r];   // k-1
                rowmap[q]=r;
            }
            else
            {
                L.p[q]=M.p[r];
                L.n[q]=M.n[r];   // i+1
                L.s[q]=M.s[r];   // i-1
                L.w[q]=M.w[r];   // j+1
                L.e[q]=M.e[r];   // j-1
                L.t[q]=M.t[r];   // k+1
                L.b[q]=M.b[r];   // k-1
            }

            const double fv=rhs.V[r], uv=f(i,j,k);
            L.f[q]=fv;
            L.u[q]=uv;
            L.act[q]=1;

            if(uv!=uv || fv!=fv)
            err=1;
        }
        else
        {
            //  solid, object or excluded cell: identity row
            if(fp32)
            {
                L.pf[q]=1.0f;
                L.nf[q]=L.sf[q]=L.wf[q]=L.ef[q]=L.tf[q]=L.bf[q]=0.0f;
            }
            else
            {
                L.p[q]=1.0;
                L.n[q]=L.s[q]=L.w[q]=L.e[q]=L.t[q]=L.b[q]=0.0;
            }

            L.f[q]=0.0;
            L.u[q]=0.0;
            L.act[q]=0;
        }
    }

    if(err)
    p->solver_error=1;

    if(fp32)
    build_colrow0(L);
}

void reefmg::fillbackvec5(lexer *p, field &f)
{
    sc_level &L=mg.fine();

    PLAINLOOP
    {
        const long q=L.idx(i,j,k);
        if(L.act[q])
        f(i,j,k)=L.u[q];
    }
}
