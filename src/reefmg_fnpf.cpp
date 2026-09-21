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

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

reefmg::reefmg(lexer *p, ghostcell *pgc, int solve_input, int precon_input)
{
    //  N 10 5x : 50 V-cycles only, 51 BiCGStab preconditioned by the V-cycle
    //  (default), 52 V-cycles with fallback to BiCGStab.
    //
    //  BiCGStab is the default because it costs the same as plain cycling on
    //  a well-shaped grid and degrades far more gracefully on an awkward one:
    //  on ragged local extents it needed 5 cycles where plain cycling needed 13.
    solve_mode = solve_input-50;
    if(solve_mode<0 || solve_mode>2) solve_mode=1;

    presweep  = 1;
    postsweep = 1;

    //  N 12 1 restores symmetric (forward and backward) line sweeps.  The
    //  default alternates one forward pass before the coarse correction and
    //  one backward pass after, which halves the smoothing cost without
    //  costing convergence on a grid-aligned anisotropy.
    mg.set_sweepstyle(p->N12==1 ? 1 : 0);

    if(precon_input>=1 && precon_input<=3)
    {
        presweep  = precon_input;
        postsweep = precon_input;
    }

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
        cout<<"semicoarsen: periodic boundaries are not supported yet - the halo "
            <<"exchange has no wraparound.  Use N 10 10-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2749);
    }

    //  The line relaxation solves each sigma column exactly, which requires
    //  the column to live on one rank.  For FNPF that is the natural layout.
    if(npz>1)
    {
        if(p->mpirank==0)
        cout<<"semicoarsen: the grid is decomposed in z ("<<npz<<" ranks).  The vertical "
            <<"line solver needs each sigma column on one rank - decompose in x and y "
            <<"only, or use N 10 10-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2750);
    }

    if(npx*npy!=nprocs)
    {
        if(p->mpirank==0)
        cout<<"semicoarsen: the decomposition is not a Cartesian product ("
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
            cout<<"semicoarsen: the decomposition has inconsistent block sizes - ranks "
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
        cout<<"semicoarsen: "<<mg.err()<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2752);
    }

    if(p->mpirank==0)
    {
        cout<<"semicoarsen: procs "<<npx<<" x "<<npy<<", local box "
            <<p->knox<<" x "<<p->knoy<<" x "<<p->knoz
            <<", multigrid levels "<<mg.levels()<<endl;

        if(p->knox%2!=0 || (p->knoy>1 && p->knoy%2!=0))
        cout<<"semicoarsen: note - odd local cell numbers give ragged coarse cells.  "
            <<"The coarse operators account for this, but even knox/knoy per rank "
            <<"still converge slightly better."<<endl;
    }
}

void reefmg::start(lexer *p, fdm *a, ghostcell *pgc, field &f, vec &rhsvec, int var)
{
    if(p->mpirank==0)
    cout<<"semicoarsen: start() not implemented - use N 10 10-19 (hypre) for this equation."<<endl;

    MPI_Abort(MPI_COMM_WORLD,-2753);
}

void reefmg::startf(lexer *p, ghostcell *pgc, field &f, vec &rhs, matrix_diag &M, int var)
{
    if(p->mpirank==0)
    cout<<"semicoarsen: startf() not implemented - use N 10 10-19 (hypre) for this equation."<<endl;

    MPI_Abort(MPI_COMM_WORLD,-2754);
}

void reefmg::startV(lexer *p, ghostcell *pgc, double *f, vec &rhs, matrix_diag &M, int var)
{
    if(p->mpirank==0)
    cout<<"semicoarsen: startV() not implemented - use N 10 10-19 (hypre) for this equation."<<endl;

    MPI_Abort(MPI_COMM_WORLD,-2755);
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
        cout<<"semicoarsen: startF() only implemented for var==8 (FNPF Laplace)."<<endl;

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
    cout<<"semicoarsen: cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres
        <<"  coarsen "<<setprecision(3)<<coarsentime
        <<" s  solve "<<setprecision(3)<<solvetime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"semicoarsen: WARNING - iteration limit N 46 = "<<p->N46<<" reached, res "
        <<relres<<endl;
}

void reefmg::fill_matrix8(lexer *p, double *f, vec &rhs, matrix_diag &M)
{
    sc_level &L=mg.fine();

    std::fill(L.p.begin(),L.p.end(),0.0);
    std::fill(L.n.begin(),L.n.end(),0.0); std::fill(L.s.begin(),L.s.end(),0.0);
    std::fill(L.w.begin(),L.w.end(),0.0); std::fill(L.e.begin(),L.e.end(),0.0);
    std::fill(L.t.begin(),L.t.end(),0.0); std::fill(L.b.begin(),L.b.end(),0.0);
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

            L.p[q]=M.p[n];
            L.n[q]=M.n[n];   // i+1
            L.s[q]=M.s[n];   // i-1
            L.w[q]=M.w[n];   // j+1
            L.e[q]=M.e[n];   // j-1
            L.t[q]=M.t[n];   // k+1
            L.b[q]=M.b[n];   // k-1

            L.f[q]=rhs.V[n];
            L.u[q]=f[FIJK];
            L.act[q]=1;

            if(L.u[q]!=L.u[q] || L.f[q]!=L.f[q])
            p->solver_error=1;
        }

        FSWDCHECK
        {
            L.p[q]=1.0;      // identity row, passes harmlessly through the line solve
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
