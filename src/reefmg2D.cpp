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

#include "reefmg2D.h"
#include "lexer.h"
#include "ghostcell.h"
#include "slice.h"
#include "matrix2D.h"
#include "vec2D.h"
#include "hypre_struct2D.h"
#include "slice4.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

reefmg2D::reefmg2D(lexer *p, ghostcell *pgc)
{
    //  Same settings and flags as the 3D connector (reefmg.cpp):
    //  BiCGStab right-preconditioned by one V-cycle, one pre- and one
    //  post-sweep, N 12 sweep style / ordering, N 14 storage precision,
    //  N 15 coarse-grid agglomeration.
    solve_mode = 1;

    presweep  = 1;
    postsweep = 1;

    mg.set_sweepstyle(p->N12==1 ? 1 : 0);
    mg.set_ordering(p->N12==2 ? 1 : 0);
    mg.set_precision(p->N14==32 ? 32 : 64);
    mg.set_agglomeration(p->N15==1 ? 1 : 0);

    Mcur=0;
    mg.set_fine_operator(this);

    topology(p,pgc);
}

reefmg2D::~reefmg2D()
{
}

void reefmg2D::topology(lexer *p, ghostcell *pgc)
{
    int nprocs;
    MPI_Comm_size(pgc->mpi_comm,&nprocs);

    int nd=0;
    int dims[3]={1,1,1}, per[3]={0,0,0}, crd[3]={0,0,0};
    MPI_Cartdim_get(pgc->cart(),&nd);
    MPI_Cart_get(pgc->cart(),nd,dims,per,crd);

    npx = dims[0];
    npy = nd>1 ? dims[1] : 1;
    const int npz = nd>2 ? dims[2] : 1;
    const int cx = crd[0];
    const int cy = nd>1 ? crd[1] : 0;

    //  one coupled horizontal direction for a 1D flume, two otherwise
    ndir = (p->gknoy>1) ? 2 : 1;

    if(p->periodic1>0 || p->periodic2>0 || p->periodic3>0)
    {
        if(p->mpirank==0)
        cout<<"REEFMG2D periodic boundaries are not supported yet - the halo "
            <<"exchange has no wraparound.  Use N 10 11-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2760);
    }

    if(npz>1)
    {
        if(p->mpirank==0)
        cout<<"REEFMG2D the grid is decomposed in z ("<<npz<<" ranks) - decompose "
            <<"in x and y only, or use N 10 11-19 (hypre)."<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2761);
    }

    //  the halo exchange needs equal knox along a process column and equal
    //  knoy along a process row, see reefmg::topology
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
            cout<<"REEFMG2D the decomposition has inconsistent block sizes - ranks "
                <<"sharing an x position must share knox, and likewise in y.  "
                <<"Use N 10 11-19 (hypre)."<<endl;

            MPI_Abort(MPI_COMM_WORLD,-2762);
        }
    }

    //  a single layer: the column solve becomes point Gauss-Seidel
    if(!mg.setup(pgc->cart(),
                 p->knox,p->knoy,1,p->gknox,p->gknoy,p->N13,
                 p->DXN+p->marge-1, p->DYN+p->marge-1))
    {
        if(p->mpirank==0)
        cout<<"REEFMG2D "<<mg.err()<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2763);
    }

    if(p->mpirank==0)
    {
        cout<<"REEFMG2D procs "<<npx<<" x "<<npy<<", local box "
            <<p->knox<<" x "<<p->knoy
            <<", multigrid levels "<<mg.levels()
            <<(p->N13==0 ? " (truncated per solve to the screening length)" : "")<<endl;

        if(mg.agglomerated())
        cout<<"REEFMG2D coarse-grid agglomeration: "<<mg.agglomerated_cells()
            <<" cells gathered on every rank"<<endl;
    }
}

//  ---------------------------------------------------------------------------
//  Levels worth keeping for a screened operator.
//
//  For  -div(a grad q) + c q  the error decays over l = sqrt(a/c).  On a level
//  with spacing H >> l the operator is diagonally dominant and Gauss-Seidel
//  alone converges quickly, so the hierarchy needs to reach H ~ l, not the
//  global grid.  ldx is the largest l/dx over all active rows; two levels
//  beyond log2(ldx) keep a margin, and the coarsest level in use still gets
//  the coarse sweeps and the Krylov wrapper, so an underestimate costs
//  iterations, never accuracy.
//  ---------------------------------------------------------------------------
int reefmg2D::screened_levels(double ldx, int nlev) const
{
    if(ldx>=1.0e20)
    return nlev;

    const int need = 2 + (int)ceil(log2(ldx>1.0 ? ldx : 1.0));

    return need<nlev ? need : nlev;
}

void reefmg2D::start(lexer *p, ghostcell *pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var)
{
    p->solveriter=0;
    starttime=pgc->timer();

    bool singular=false;
    const double ldx=fill(p,pgc,mg,f,M,rhsvec,singular);

    //  pure Neumann: the potential-flow initialisation (sflow_potential_f)
    if(singular)
    {
        solve_neumann(p,pgc,f,M,xvec,rhsvec,var);
        return;
    }

    if(p->N13==0)
    mg.set_active_levels(screened_levels(ldx,mg.levels()));

    mg.coarsen();

    double relres=1.0;
    p->solveriter = mg.solve_auto(p->N44,p->N46,relres,presweep,postsweep,solve_mode);
    p->final_res=relres;

    //  A right-hand side that overflows the norm (a flow that has already blown
    //  up) ends the iteration at once with the initial guess, as hypre does;
    //  NaN in the right-hand side or the initial guess is reported through
    //  solver_error by fill().

    fillback(p,mg,f);

    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"REEFMG2D levels "<<mg.active_levels()<<"/"<<mg.levels()
        <<"  l/dx "<<setprecision(3)<<(ldx<1.0e20 ? ldx : -1.0)
        <<"  cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres
        <<"  "<<setprecision(3)<<pgc->timer()-starttime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"REEFMG2D WARNING - iteration limit N 46 = "<<p->N46<<" reached, res "
        <<relres<<endl;
}

//  ---------------------------------------------------------------------------
//  matrix2D -> fine level
//
//  Rows are numbered in SLICELOOP4 order, as in every assembly that calls a
//  solver2D (sflow_pjm_lin/quad::poisson, sflow_potential_f::laplace) and as
//  in hypre_struct2D::fill_matrix.  Stencil: n i+1, s i-1, w j+1, e j-1.
//
//  Three kinds of rows:
//   active     : at least one off-diagonal - copied into the multigrid.
//   decoupled  : no off-diagonals - hydrostatic, breaking or dry cells, which
//                the pressure assembly turns into q = 0 rows.  Solved here
//                as q = rhs/p and kept out of the multigrid.
//   excluded   : zero diagonal, or no matrix row (flagslice4 < 0).
//  A coupling from an active row to anything that is not active is moved to
//  the right-hand side with that cell's value and dropped from the stencil,
//  as the CFD assembly does for walls.  The diagonal keeps it, so the coarse
//  operators - which average the row-sum excess over the agglomerate - see a
//  Dirichlet boundary there instead of a neighbour.  The same drop is
//  applied in fine_apply(), so Krylov iteration and preconditioner agree.
//  ---------------------------------------------------------------------------
double reefmg2D::fill(lexer *p, ghostcell *pgc, reefmg_core &core, slice &f, matrix2D &M,
                      vec2D &rhs, bool &singular)
{
    sc_level &L=core.fine();
    const bool fp32=(core.precision()==32);
    const long N=L.size();
    const long sx=(long)(L.ny+2), sy=1;

    Mcur=&M;

    //  row numbering - redone every solve, flagslice4 may change with the bed
    cval.assign((size_t)p->imax*(size_t)p->jmax,-1);
    int cnt=0;
    SLICELOOP4
    {
        cval[IJ]=cnt;
        ++cnt;
    }

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

    rowmap.assign(N,-1);
    nbm.assign(N,0);
    dval.assign(N,0.0);
    dflag.assign(N,0);
    amask.assign(N,0.0);
    dhalo.assign(N,0.0);

    int err=0;

    //  pass 1: classify the rows
    ILOOP
    JLOOP
    {
        if(p->flagslice4[IJ]<=0)
        continue;

        const long q=L.idx(i,j,0);
        const int r=cval[IJ];
        const double mp=M.p[r];

        if(mp==0.0)
        continue;

        if(M.n[r]==0.0 && M.s[r]==0.0 && M.w[r]==0.0 && M.e[r]==0.0)
        {
            const double v=rhs.V[r]/mp;
            dflag[q]=1;
            dval[q]=dhalo[q]=v;

            if(v!=v)
            err=1;

            continue;
        }

        rowmap[q]=r;
        amask[q]=1.0;
    }

    //  what the neighbours on the other ranks are; the global edge stays 0,
    //  i.e. outside the problem
    core.halo_vec(L,amask);
    core.halo_vec(L,dhalo);

    //  pass 2: copy the active rows
    double ldx=0.0, excmax=0.0;

    ILOOP
    JLOOP
    {
        const long q=L.idx(i,j,0);
        const int r=rowmap[q];

        if(r<0)
        {
            //  identity row, passes harmlessly through the smoother
            if(fp32) L.pf[q]=1.0f;
            else     L.p[q]=1.0;
            continue;
        }

        const double mp=M.p[r];
        double c[4]={M.n[r],M.s[r],M.w[r],M.e[r]};
        const long nq[4]={q+sx,q-sx,q+sy,q-sy};

        //  screening from the assembled row: reaction = row-sum excess,
        //  diffusion = sum of the couplings
        const double react=mp+c[0]+c[1]+c[2]+c[3];
        const double diff=fabs(c[0])+fabs(c[1])+fabs(c[2])+fabs(c[3]);

        if(react>1.0e-12*fabs(mp))
        {
            const double l=sqrt(diff/(2.0*ndir*react));
            if(l>ldx) ldx=l;
        }
        else
        ldx=1.0e30;

        double fv=rhs.V[r];
        unsigned char m=0;

        for(int d=0;d<4;++d)
        {
            if(c[d]==0.0)
            continue;

            if(amask[nq[d]]>0.5)
            m|=(unsigned char)(1<<d);
            else
            {
                fv-=c[d]*dhalo[nq[d]];
                c[d]=0.0;
            }
        }

        //  row-sum excess of what the multigrid sees: reaction plus any
        //  Dirichlet coupling just dropped - zero everywhere means pure Neumann
        const double exc=(mp+c[0]+c[1]+c[2]+c[3])/fabs(mp);
        if(exc>excmax) excmax=exc;

        if(fp32)
        {
            L.pf[q]=(float)mp;
            L.nf[q]=(float)c[0];   // i+1
            L.sf[q]=(float)c[1];   // i-1
            L.wf[q]=(float)c[2];   // j+1
            L.ef[q]=(float)c[3];   // j-1
        }
        else
        {
            L.p[q]=mp;
            L.n[q]=c[0];
            L.s[q]=c[1];
            L.w[q]=c[2];
            L.e[q]=c[3];
        }

        const double uv=f(i,j);
        L.f[q]=fv;
        L.u[q]=uv;
        L.act[q]=1;
        nbm[q]=m;

        if(uv!=uv || fv!=fv)
        err=1;
    }

    if(err)
    p->solver_error=1;

    excmax=pgc->globalmax(excmax);
    singular=(excmax<1.0e-10);

    return pgc->globalmax(ldx);
}

void reefmg2D::fillback(lexer *p, reefmg_core &core, slice &f)
{
    sc_level &L=core.fine();

    ILOOP
    JLOOP
    {
        const long q=L.idx(i,j,0);

        if(L.act[q])
        f(i,j)=L.u[q];

        else if(dflag[q])
        f(i,j)=dval[q];
    }
}

//  y = A x with the exact double operator, read from matrix2D, with the same
//  couplings dropped as in fill().  Used by the Krylov iteration in fp32 mode.
void reefmg2D::fine_apply(const sc_level &L,const double *x,double *y)
{
    const matrix2D &M=*Mcur;
    const long sx=(long)(L.ny+2), sy=1;

    for(int ii=0;ii<L.nx;++ii)
    for(int jj=0;jj<L.ny;++jj)
    {
        const long q=L.idx(ii,jj,0);
        const int r=rowmap[q];

        if(r<0)
        {
            y[q]=x[q];
            continue;
        }

        const unsigned char m=nbm[q];
        double v=M.p[r]*x[q];

        if(m&1) v+=M.n[r]*x[q+sx];
        if(m&2) v+=M.s[r]*x[q-sx];
        if(m&4) v+=M.w[r]*x[q+sy];
        if(m&8) v+=M.e[r]*x[q-sy];

        y[q]=v;
    }
}

//  hypre GMRES preconditioned by PFMG, as with N 10 12 / N 11 11
void reefmg2D::solve_hypre(lexer *p, ghostcell *pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var)
{
    const int n10=p->N10, n11=p->N11;
    p->N10=12;
    p->N11=11;
    {
        hypre_struct2D fallback(p,pgc);
        fallback.start(p,pgc,f,M,xvec,rhsvec,var);
    }
    p->N10=n10;
    p->N11=n11;
}

//  ---------------------------------------------------------------------------
//  Pure-Neumann problems: sflow_potential_f, the depth-averaged Laplace for
//  the discharge potential with walls, dry cells and prescribed in- and
//  outflow fluxes.  As for the NHFLOW potential (reefmg::start_solver44), the
//  smooth mode is unscreened and needs the whole hierarchy and a coarsest
//  level that resolves it, so this runs once on its own temporary double-
//  precision, full-depth, agglomerated hierarchy.  If that cannot be had, or
//  the solve diverges, hypre GMRES+PFMG solves it.
//  ---------------------------------------------------------------------------
static const long POT2D_COARSEST_MAX  = 512;
static const int  POT2D_COARSE_SWEEPS = 128;

void reefmg2D::solve_neumann(lexer *p, ghostcell *pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var)
{
    reefmg_core pot;
    pot.set_precision(64);
    pot.set_coarse_sweeps(POT2D_COARSE_SWEEPS);
    pot.set_agglomeration(1);

    if(!pot.setup(pgc->cart(),
                  p->knox,p->knoy,1,p->gknox,p->gknoy,0,
                  p->DXN+p->marge-1, p->DYN+p->marge-1))
    {
        if(p->mpirank==0)
        cout<<"REEFMG2D pure-Neumann solve: "<<pot.err()<<endl;

        MPI_Abort(MPI_COMM_WORLD,-2764);
    }

    const sc_level &C=pot.coarsest();

    if(!pot.agglomerated() && (long)C.gnx*C.gny > POT2D_COARSEST_MAX)
    {
        if(p->mpirank==0)
        cout<<"REEFMG2D pure-Neumann solve: coarsest grid "<<C.gnx<<" x "<<C.gny
            <<" too large to resolve reliably and to agglomerate - solving it with "
            <<"hypre GMRES+PFMG instead."<<endl;

        solve_hypre(p,pgc,f,M,xvec,rhsvec,var);
        return;
    }

    //  initial guess, kept in case the solve has to be handed over
    slice4 f0(p);
    SLICELOOP4
    f0(i,j)=f(i,j);

    bool singular=true;
    fill(p,pgc,pot,f,M,rhsvec,singular);
    pot.coarsen();

    double relres=1.0;
    p->solveriter=pot.solve_auto(p->N44,p->N46,relres,1,1,1);
    p->final_res=relres;

    //  An inconsistent pure-Neumann system - prescribed in- and outflow that
    //  do not balance - has no solution, and BiCGStab then drifts off instead
    //  of stopping.  Do not hand that back: restart from the initial guess
    //  with hypre, i.e. behave exactly as without REEFMG.
    if(!(relres<=1.0))
    {
        if(p->mpirank==0)
        cout<<"REEFMG2D pure-Neumann solve diverged (res "<<relres<<") - the system is "
            <<"probably inconsistent (in- and outflow do not balance).  Handing it to "
            <<"hypre GMRES+PFMG."<<endl;

        SLICELOOP4
        f(i,j)=f0(i,j);

        solve_hypre(p,pgc,f,M,xvec,rhsvec,var);
        return;
    }

    fillback(p,pot,f);

    if(p->mpirank==0)
    cout<<"REEFMG2D pure-Neumann solve: levels "<<pot.levels()
        <<(pot.agglomerated()? " + agglomerated coarse grid" : "")
        <<"  cycles "<<p->solveriter
        <<"  res "<<setprecision(3)<<relres<<"  "<<setprecision(3)
        <<pgc->timer()-starttime<<" s"<<endl;

    if(p->solveriter>=p->N46 && p->mpirank==0)
    cout<<"REEFMG2D pure-Neumann solve: WARNING - iteration limit N 46 reached.  "
        <<"If the prescribed inflow and outflow do not balance, the system is "
        <<"inconsistent and cannot converge with any solver."<<endl;
}
