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

#include"nhflow_amr.h"
#include"nhflow_amr_fill.h"
#include"nhflow_momentum_SI.h"
#include"nhflow_longrun.h"
#include"si_eta2D.h"
#include"reefmg2D.h"
#include"vec2D.h"
#include"matrix2D.h"
#include"reefamr_krylov.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"nhflow_curv.h"
#include"sflow_curv.h"
#include<mpi.h>
#include<iomanip>

// Hydrostatic semi-implicit free surface (A 534 1) on the multi-block grid (G 50 1, DIVEMesh R 20 1).
//
// Every block (the native grid -1 and the blocks 1.. of lev[0]) has its own nhflow_momentum_SI. The
// step runs its parts on all blocks; between them the cells around the blocks and the ring of the
// native grid take the values of the cells they stand for, as the stages of the explicit integrators
// (fill_level0, fill_run0):
//   si_restore                       the water column of the step before
//   fill_level0, si_fillF            eta, WL, wet, U, V, W, UH, VH, WH, omega of t^n around the blocks
//   si_begin                         t^n, face depths, bed friction, open boundaries
//   per predictor stage:  si_stage_vel, fill U, V, W, si_stage_rhs
//   si_predictor_end, fill UHD, VHD, si_assemble   the eta system of every block (si_eta2D::build)
//   si_solve                         one eta system over the own cells of all blocks: BiCGStab
//                                    (reefamr_bicgstab), preconditioned by one REEFMG2D V-cycle per
//                                    grid on its own rows (restricted additive Schwarz); the cells
//                                    around the blocks and the ring are filled from the vectors in
//                                    every product (xslice), so a face between two blocks couples
//                                    them as an interior face; the wetting and drying iteration over
//                                    all blocks together
//   fill eta^{n+1}, si_end           faces, WL, cells, omega
// A face between two blocks is computed by both from the same values (the filled cells hold the
// values of the cells they stand for), its discharge is the same on both sides.

using namespace nhflow_amr_detail;

namespace
{
typedef reefamr_comms_off comms_off;

// the vector space of the composite eta solve for reefamr_bicgstab
struct si_space
{
    nhflow_amr *a;
    static const int B=0, R=1, RH=2, PV=3, VV=4, S=5, T=6, PH=7, SH=8;
    void apply(int x, int y) { a->si_apply(x,y); }
    void prec(int r, int z) { a->si_prec(r,z); }
    double dot(int x, int y) { return a->si_dot(x,y); }
    void op_start() { a->si_opstart(); }
    void op_p(double beta, double om) { a->si_op_p(beta,om); }
    void op_s(double alp) { a->si_op_s(alp); }
    void op_x(double alp, double om) { a->si_op_x(alp,om); }
    void op_restart() { a->si_restart(); }
};

const int SIX = 9;      // the unknown (eta^{n+1})
const int SINV = 10;
}

nhflow_momentum_SI* nhflow_amr::gsi(int g)
{
    return static_cast<nhflow_momentum_SI*>(gmom(g));
}

// a cell of grid g that is an unknown of the composite system: a fluid cell of the native grid outside
// the ring, a fluid cell of a block inside its box
bool nhflow_amr::si_own(int g, int ii, int jj)
{
    lexer *q = glex(g);

    if(q->flagslice4[lij(q,ii,jj)]<0)
    return false;

    if(g<0)
    return !ring0(ii,jj);

    reefamr_patch *c = P[g];
    return ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny;
}

void nhflow_amr::si_setup()
{
    sigrid.clear();
    sigrid.push_back(-1);

    for(int id : lev[0])
    sigrid.push_back(id);

    sirow.assign(sigrid.size(),vector<int>());
    sivec.assign(sigrid.size(),vector<slice*>());

    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        const int g = sigrid[gi];
        lexer *q = glex(g);

        // the rows in the order of SLICELOOP4 (si_eta2D::build)
        int n=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(q->flagslice4[lij(q,ii,jj)]>0)
        {
            if(si_own(g,ii,jj))
            {
            sirow[gi].push_back(ii);
            sirow[gi].push_back(jj);
            sirow[gi].push_back(n);
            }
        ++n;
        }

        sinrow.push_back(n);

        for(int v=0; v<SINV; ++v)
        {
        slice *s = new slice4(q);
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        (*s)(ii,jj) = 0.0;
        sivec[gi].push_back(s);
        }
    }

    // P7: the V-cycle preconditioner of every grid (native: over its ranks; blocks: on their rank)
    simg.assign(sigrid.size(),nullptr);
    sir.assign(sigrid.size(),nullptr);
    siz.assign(sigrid.size(),nullptr);
    sirebuild.assign(sigrid.size(),1);

    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        const int g = sigrid[gi];
        lexer *q = glex(g);

        if(g<0)
        simg[gi] = new reefmg2D(q,sipgc);

        else
        {
            comms_off guard(sipgc);
            pscope ps(sipgc,NP(g)->d,d0);
            simg[gi] = new reefmg2D(q,sipgc,true);
        }

        sir[gi] = new vec2D(q,max(sinrow[gi],1));
        siz[gi] = new slice4(q);
    }
}

// f(g, lexer, fdm, SI object, stage objects) on every block: the native grid with the MPI exchanges of
// its partition, the blocks 1.. on their rank alone (pscope)
template<class F>
void nhflow_amr::si_each(ghostcell *pgc, F f)
{
    f(-1,p0,d0,gsi(-1),*S0p);

    comms_off guard(pgc);
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        pscope ps(pgc,c->d,d0);
        f(id,c->pp,c->d,gsi(id),c->S);
    }
}

// nf layer arrays A(g,m) of every block: the cells around the blocks and the ring of the native grid
// from the cells they stand for, then the partition halo of the native grid
template<class A>
void nhflow_amr::si_fill3(ghostcell *pgc, int nf, int tag, A arr)
{
    const int K = p0->knoz;
    const int nv = nf*K;

    fill_run0(nv,tag,
              [&](const reefamr_fill &f, double *v)
              {
                  if(f.kind!=0)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = 0.0;
                      return;
                  }
                  lexer *q = glex(f.g);
                  for(int m=0; m<nf; ++m)
                  {
                      const double *src = arr(f.g,m);
                      for(int k=0; k<K; ++k)
                      v[m*K+k] = src[cidx(q,f.si,f.sj,k)];
                  }
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  lexer *pp = q->pp;
                  for(int m=0; m<nf; ++m)
                  {
                      double *dst = arr(id,m);
                      for(int k=0; k<K; ++k)
                      dst[cidx(pp,f.di,f.dj,k)] = w[m*K+k];
                  }
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  for(int m=0; m<nf; ++m)
                  {
                      double *dst = arr(-1,m);
                      for(int k=0; k<K; ++k)
                      dst[cidx(p0,f.di,f.dj,k)] = w[m*K+k];
                  }
              });

    if(p0->mpi_size>1)
    for(int m=0; m<nf; ++m)
    {
        pgc->gcparaxV(p0,arr(-1,m),14);
        pgc->gcparacoxV(p0,arr(-1,m),14);
    }
}

// omega (sigma velocity, K+1 levels) of all blocks: the cells around the blocks and the ring of the native
// grid from the cells they stand for (the face columns between two blocks take the omega of both cells,
// col_coef), then the partition halo of the native grid
void nhflow_amr::si_fillF(ghostcell *pgc)
{
    const int K = p0->knoz;
    const int nv = K+1;

    fill_run0(nv,7640,
              [&](const reefamr_fill &f, double *v)
              {
                  if(f.kind!=0)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = 0.0;
                      return;
                  }
                  lexer *q = glex(f.g);
                  const double *src = gfd(f.g)->omegaF;
                  for(int k=0; k<=K; ++k)
                  v[k] = src[fidx(q,f.si,f.sj,k)];
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
              {
                  lexer *pp = q->pp;
                  double *dst = NP(q)->d->omegaF;
                  for(int k=0; k<=K; ++k)
                  dst[fidx(pp,f.di,f.dj,k)] = w[k];
              },
              [&](const reefamr_fill &f, const double *w)
              {
                  for(int k=0; k<=K; ++k)
                  d0->omegaF[fidx(p0,f.di,f.dj,k)] = w[k];
              });

    if(p0->mpi_size>1)
    pgc->start7S(p0,d0->omegaF,17);
}

// Krylov vector v of all blocks: the cells around the blocks, the ring and the halo of the native grid
void nhflow_amr::si_fillv(ghostcell *pgc, int v)
{
    xslice(1,7620,[&](int g, int k) -> slice*
    {
        for(size_t gi=0; gi<sigrid.size(); ++gi)
        if(sigrid[gi]==g)
        return sivec[gi][v];
        return sivec[0][v];
    });
}

// the water level boundaries of every block: the boundary value in the ghost cell of eta^{n+1}
void nhflow_amr::si_dirichlet(ghostcell *pgc)
{
    for(int g : sigrid)
    {
        nhflow_momentum_SI *m = gsi(g);

        for(const auto &f : m->bcf)
        if(f.type==nhflow_momentum_SI::bc_level)
        m->etaI(f.gi,f.gj) = m->etaB(f.gi,f.gj);
    }
}

void nhflow_amr::si_apply(int x, int y)
{
    si_fillv(sipgc,x);

    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        const matrix2D &M = gsi(sigrid[gi])->psi->mat();
        slice &X = *sivec[gi][x];
        slice &Y = *sivec[gi][y];
        const vector<int> &r = sirow[gi];

        for(size_t q=0; q<r.size(); q+=3)
        {
            const int ii=r[q], jj=r[q+1], n=r[q+2];

            Y(ii,jj) = M.p[n]*X(ii,jj) + M.n[n]*X(ii+1,jj) + M.s[n]*X(ii-1,jj)
                     + M.w[n]*X(ii,jj+1) + M.e[n]*X(ii,jj-1);
        }
    }
}

void nhflow_amr::si_prec(int r, int z)
{
    // one V-cycle per grid on the residual of its own rows (zero on the cells around it), the result
    // on its own rows: restricted additive Schwarz over the grids with a multigrid subdomain solve
    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        const int g = sigrid[gi];
        lexer *q = glex(g);
        matrix2D &M = const_cast<matrix2D&>(gsi(g)->psi->mat());
        slice &R = *sivec[gi][r];
        slice &Z = *sivec[gi][z];
        vec2D &rv = *sir[gi];
        const vector<int> &w = sirow[gi];

        std::fill(rv.V.begin(),rv.V.end(),0.0);
        for(size_t k=0; k<w.size(); k+=3)
        rv.V[w[k+2]] = R(w[k],w[k+1]);

        if(g<0)
        simg[gi]->vcycle(q,sipgc,*siz[gi],M,rv,sirebuild[gi]);

        else
        {
            comms_off guard(sipgc);
            pscope ps(sipgc,NP(g)->d,d0);
            simg[gi]->vcycle(q,sipgc,*siz[gi],M,rv,sirebuild[gi]);
        }

        sirebuild[gi] = 0;

        for(size_t k=0; k<w.size(); k+=3)
        Z(w[k],w[k+1]) = (*siz[gi])(w[k],w[k+1]);
    }
}

double nhflow_amr::si_dot(int a, int b)
{
    double s=0.0;

    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        slice &A = *sivec[gi][a];
        slice &B = *sivec[gi][b];
        const vector<int> &w = sirow[gi];

        for(size_t q=0; q<w.size(); q+=3)
        s += A(w[q],w[q+1])*B(w[q],w[q+1]);
    }

    double r=0.0;
    MPI_Allreduce(&s,&r,1,MPI_DOUBLE,MPI_SUM,sipgc->mpi_comm);
    return r;
}

// the vector operations of reefamr_bicgstab over the own cells
namespace
{
template<class F>
void si_rows(vector<vector<slice*>> &V, vector<vector<int>> &rows, F f)
{
    for(size_t gi=0; gi<V.size(); ++gi)
    {
        const vector<int> &w = rows[gi];
        for(size_t q=0; q<w.size(); q+=3)
        f(V[gi],w[q],w[q+1]);
    }
}
}

void nhflow_amr::si_opstart()
{
    si_rows(sivec,sirow,[](vector<slice*> &v, int ii, int jj)
    {
        const double r = (*v[si_space::B])(ii,jj) - (*v[si_space::VV])(ii,jj);
        (*v[si_space::R])(ii,jj) = r;
        (*v[si_space::RH])(ii,jj) = r;
        (*v[si_space::PV])(ii,jj) = 0.0;
        (*v[si_space::VV])(ii,jj) = 0.0;
    });
}

void nhflow_amr::si_op_p(double beta, double om)
{
    si_rows(sivec,sirow,[&](vector<slice*> &v, int ii, int jj)
    {
        (*v[si_space::PV])(ii,jj) = (*v[si_space::R])(ii,jj) + beta*((*v[si_space::PV])(ii,jj) - om*(*v[si_space::VV])(ii,jj));
    });
}

void nhflow_amr::si_op_s(double alp)
{
    si_rows(sivec,sirow,[&](vector<slice*> &v, int ii, int jj)
    {
        (*v[si_space::S])(ii,jj) = (*v[si_space::R])(ii,jj) - alp*(*v[si_space::VV])(ii,jj);
    });
}

void nhflow_amr::si_op_x(double alp, double om)
{
    si_rows(sivec,sirow,[&](vector<slice*> &v, int ii, int jj)
    {
        (*v[SIX])(ii,jj) += alp*(*v[si_space::PH])(ii,jj) + om*(*v[si_space::SH])(ii,jj);
        (*v[si_space::R])(ii,jj) = (*v[si_space::S])(ii,jj) - om*(*v[si_space::T])(ii,jj);
    });
}

void nhflow_amr::si_restart()
{
    si_rows(sivec,sirow,[](vector<slice*> &v, int ii, int jj)
    {
        (*v[si_space::RH])(ii,jj) = (*v[si_space::R])(ii,jj);
        (*v[si_space::PV])(ii,jj) = 0.0;
        (*v[si_space::VV])(ii,jj) = 0.0;
    });
}

// the eta system of all blocks with the wetting and drying iteration of si_eta2D::solve_wetdry
void nhflow_amr::si_solve(lexer *p, ghostcell *pgc)
{
    const double t0 = pgc->timer();
    const int itmax = p->A537==1 ? 20 : 1;
    int it_done=0, iter=0;

    for(int g : sigrid)
    gsi(g)->psi->wd_start(glex(g));

    si_space sp{this};

    for(int it=0; it<itmax; ++it)
    {
        // the rows of every block, the right-hand side and the initial guess
        for(size_t gi=0; gi<sigrid.size(); ++gi)
        {
            const int g = sigrid[gi];
            lexer *q = glex(g);
            nhflow_momentum_SI *m = gsi(g);
            si_eta2D *e = m->psi;

            e->wd_rows(q,m->rhsW,gfd(g)->depth,m->bdiag);
            e->build(q,m->ax,m->ay,e->wd_rhs(),m->etaI,e->wd_react());

            const vec2D &b = e->rhsv();
            const vector<int> &w = sirow[gi];
            for(size_t k=0; k<w.size(); k+=3)
            {
            (*sivec[gi][si_space::B])(w[k],w[k+1]) = b.V[w[k+2]];
            (*sivec[gi][SIX])(w[k],w[k+1]) = m->etaI(w[k],w[k+1]);
            }
        }

        // the matrices changed: the preconditioner coarsens again at its first use (keeping the
        // hierarchy over the passes of the wetting and drying iteration doubles the iterations)
        std::fill(sirebuild.begin(),sirebuild.end(),1);

        si_apply(SIX,si_space::VV);

        double bn,rn;
        int restarts=0;
        iter = reefamr_bicgstab(sp,p->N44,p->N46,bn,rn,nullptr,nullptr,&restarts);

        if(iter>=p->N46 && p->mpirank==0)
        cout<<"NHFLOW A 534 1 multi-block: eta solve at the iteration limit N 46 = "<<p->N46<<", residual "<<(bn>0.0?rn/bn:rn)<<endl;

        // eta^{n+1} of the own cells, then around the blocks and in the ring
        for(size_t gi=0; gi<sigrid.size(); ++gi)
        {
            nhflow_momentum_SI *m = gsi(sigrid[gi]);
            const vector<int> &w = sirow[gi];
            for(size_t k=0; k<w.size(); k+=3)
            m->etaI(w[k],w[k+1]) = (*sivec[gi][SIX])(w[k],w[k+1]);
        }

        xslice(1,7630,[&](int g, int k) -> slice* { return &gsi(g)->etaI; });
        si_dirichlet(pgc);

        ++it_done;

        if(itmax==1)
        break;

        // new wet pattern
        int change=0;
        for(int g : sigrid)
        {
            lexer *q = glex(g);
            nhflow_momentum_SI *m = gsi(g);
            change += m->psi->wd_update(q,m->ax,m->ay,m->bdiag,m->etaI,gfd(g)->depth);
        }

        int all=0;
        MPI_Allreduce(&change,&all,1,MPI_INT,MPI_SUM,pgc->mpi_comm);

        if(all==0)
        break;
    }

    // dry cells
    int dry=0, alldry=0;
    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        si_eta2D *e = gsi(sigrid[gi])->psi;
        const vector<int> &w = sirow[gi];
        for(size_t k=0; k<w.size(); k+=3)
        if(e->wd_pattern(w[k],w[k+1])<0.5)
        ++dry;
    }
    MPI_Allreduce(&dry,&alldry,1,MPI_INT,MPI_SUM,pgc->mpi_comm);

    const double t1 = pgc->timer()-t0;

    for(int g : sigrid)
    {
        nhflow_momentum_SI *m = gsi(g);
        m->psi_iter = iter;
        m->psi_time = t1;
        m->wd_iter = it_done;
        m->wd_cells = alldry;
    }
}

void nhflow_amr::si_step(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_momentum_func *mom, nhflow_stage_obj &S0)
{
    S0p = &S0;
    sipgc = pgc;

    if(sigrid.empty())
    si_setup();

    // the time of the blocks
    for(int id : lev[0])
    {
        lexer *pp = NP(id)->pp;
        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;
        pp->amr_follow(*p);
    }

    // the water column of the step before, then the state of t^n around the blocks
    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_restore(q,dd,pgc);
    });

    fill_level0(pgc,-1);
    si_fillF(pgc);

    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_begin(q,dd,pgc,S);
    });

    // predictor: the velocities of a stage around the blocks before its fluxes
    const int ns = gsi(-1)->si_stages();

    for(int s=0; s<ns; ++s)
    {
        si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
        {
            m->si_stage_vel(q,dd,pgc,S,s);
        });

        si_fill3(pgc,3,7601+s,[&](int g, int c) -> double*
        {
            fdm_nhf *dd = gfd(g);
            return c==0 ? dd->U : (c==1 ? dd->V : dd->W);
        });

        si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
        {
            m->si_stage_rhs(q,dd,pgc,S,s);
        });
    }

    // the cell predictor around the blocks, then the eta system of every block
    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_predictor_end(q,dd,pgc);
    });

    si_fill3(pgc,2,7610,[&](int g, int c) -> double*
    {
        return c==0 ? gsi(g)->UHD : gsi(g)->VHD;
    });

    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_assemble(q,dd,pgc);
    });

    // one eta system over all blocks
    si_solve(p,pgc);

    // faces, water level, cells
    si_each(pgc,[&](int g, lexer *q, fdm_nhf *dd, nhflow_momentum_SI *m, nhflow_stage_obj &S)
    {
        m->si_end(q,dd,pgc,S);
    });

    // the state of t^{n+1} around the blocks and in the ring (output, time step)
    fill_level0(pgc,-1);
    si_fillF(pgc);

    // the log of all blocks: volume of the own cells, discharges over the open boundaries
    double loc[4] = {0.0,0.0,0.0,0.0}, all[4];

    for(size_t gi=0; gi<sigrid.size(); ++gi)
    {
        const int g = sigrid[gi];
        lexer *q = glex(g);
        fdm_nhf *dd = gfd(g);
        nhflow_momentum_SI *m = gsi(g);
        const vector<int> &w = sirow[gi];

        for(size_t k=0; k<w.size(); k+=3)
        loc[0] += dd->WL(w[k],w[k+1])*m->carea(q,w[k],w[k+1]);

        // the native grid has summed its discharges over the ranks already
        if(g>=0)
        {
        loc[1] += m->Qin;
        loc[2] += m->Qout;
        loc[3] += m->vrelax;
        }
    }

    MPI_Allreduce(loc,all,4,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    nhflow_momentum_SI *m0 = gsi(-1);
    m0->write_log(p,all[0],all[1]+m0->Qin,all[2]+m0->Qout,all[3]+m0->vrelax,m0->psi_iter,m0->wd_iter,m0->wd_cells);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"SI eta (multi-block): iter "<<m0->psi_iter<<"  time "<<setprecision(3)<<m0->psi_time;

    if(p->A537==1)
    cout<<"  wetdry iterations "<<m0->wd_iter<<"  dry cells "<<m0->wd_cells;

    cout<<endl;
    }
}

// the time step of the semi-implicit step over the block cells (nhflow_timestep_SI::start): the
// advective limit per cell and the explicit horizontal diffusion of the predictor
void nhflow_amr::dt_si_cmin(double cfl, double &cmin, double &hmax)
{
    for(auto q : P)
    {
        nhflow_amr_patch *c = NP(q);
        lexer *pp = c->pp;
        fdm_nhf *dd = c->d;
        double *EVh = (pp->A561>0) ? dd->EVH : dd->EV;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        {
            const int n2 = lij(pp,ii,jj);

            if(pp->flagslice4[n2]<0 || pp->wet[n2]!=1 || dd->WL(ii,jj)<=pp->A544)
            continue;

            hmax = MAX(hmax,dd->WL(ii,jj));

            for(int kk=0; kk<pp->knoz; ++kk)
            {
                const int n = cidx(pp,ii,jj,kk);
                double sigma;

                if(pp->geo_curv==1)
                {
                const double sp = sqrt(dd->U[n]*dd->U[n] + dd->V[n]*dd->V[n]);
                sigma = sp/nhflow_curv::width1(pp,ii,jj) + sp/nhflow_curv::width2(pp,ii,jj);
                }

                else
                {
                sigma = fabs(dd->U[n])/pp->DXN[ii+marge];

                if(pp->j_dir==1)
                sigma += fabs(dd->V[n])/pp->DYN[jj+marge];
                }

                if(sigma>1.0e-20)
                cmin = MIN(cmin,1.0/sigma);

                if(pp->A512>0 || pp->A560>0)
                {
                const double visc = dd->VISC[n] + EVh[n];
                const double lam = pp->geo_curv==1 ? sflow_curv::inv2(pp,ii,jj)
                                 : 1.0/(pp->DXN[ii+marge]*pp->DXN[ii+marge]) + pp->y_dir/(pp->DYN[jj+marge]*pp->DYN[jj+marge]);

                cmin = MIN(cmin, 0.25/(visc*lam + 1.0e-20)/(2.0*cfl));
                }
            }
        }
    }
}

// long-run statistics (P 160): the blocks 1.. (their own cells)
void nhflow_amr::longrun_grids(nhflow_longrun *lr)
{
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        reefamr_patch *q = P[id];
        lr->add_grid(c->pp,c->d,EXT,EXT+q->nx,EXT,EXT+q->ny);
    }
}
