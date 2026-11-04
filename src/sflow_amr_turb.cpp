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

#include "fdm2D.h"
#include "ghostcell.h"
#include "increment.h"
#include "lexer.h"
#include "matrix2D.h"
#include "sflow_amr.h"
#include "sflow_bicgstab.h"
#include "sflow_turb_io.h"
#include "sflow_turb_ke_IM1.h"
#include "sflow_turb_kw_IM1.h"
#include "solver2D.h"
#include "vec2D.h"
#include <climits>
#include <mpi.h>

// Turbulence k-epsilon (A 260 1) and k-omega (A 260 2) on the patches of the mesh refinement (Cartesian grids
// without multi-block; as NHFLOW, nhflow_amr_turb.cpp).
//
// Every patch has its own model (sflow_turb_ke_IM1 / sflow_turb_kw_IM1 on the lexer and the fdm2D of the patch): its
// implicit convection and diffusion, its sources and its BiCGStab.  After the step of level 0 (sflow_f::mainloop: the
// momentum of all grids, then the model of level 0):
//   1. the face velocities of the patches from their end state (the production, the convection)
//   2. levels 1..: the cells around the patches (k, epsilon/omega, the eddy viscosity) from a patch of the same level
//      or linear from the parent (minmod slopes: k and epsilon/omega stay positive), then the patches step.  The
//      solver of a patch takes the cells around it as Dirichlet values (turb_solver: their coefficients move to the
//      right-hand side; the walls as on level 0)
//   3. the covered cells take the mean of their children, finest first; the ghost cells of level 0
//   4. the old time level of the patches (level 0: the main loop)
// The patches step with the time step of level 0, also with G 7 1 (one turbulence step per level-0 step: the eddy
// viscosity of the patches is that of the end of the last level-0 step during their steps, as on level 0).
// A fresh patch of a regrid takes k, epsilon/omega and the eddy viscosity of the old patch of the level where there
// was one, else linear from its parent (turb_regrid).
//
// A patch cut at a rank box (or placed, G 40) solves its pieces apart, with the cells of the other pieces at the start
// of the step: the turbulence fields differ from those of one rank, as an implicit solve on several ranks does at the
// solver tolerance on level 0, but by more (the lag of one step at the cuts).

namespace {
typedef reefamr_comms_off comms_off;

inline int lij(const lexer *q, int ii, int jj) { return (ii - q->imin) * q->jmax + (jj - q->jmin); }

// a patch kernel: no exchange, local reductions, no screen output of the model
struct tscope
{
    ghostcell *g;
    lexer *q;
    bool oc, ol;
    int p12, cnt;
    tscope(ghostcell *gg, lexer *pp) : g(gg), q(pp), p12(pp->P12), cnt(pp->count)
    {
        oc = g->set_comms(false);
        ol = g->set_local(true);
        q->P12 = INT_MAX;
        q->count = 1;
    }
    ~tscope()
    {
        q->P12 = p12;
        q->count = cnt;
        g->set_local(ol);
        g->set_comms(oc);
    }
};

// the solver of a patch: the coefficients of the cells beyond the computed cells that a grid fills (fluid cells)
// move to the right-hand side with their values; walls as on level 0
class sflow_amr_tsolver : public solver2D, public increment
{
public:
    sflow_amr_tsolver(lexer *p, ghostcell *pgc) : bs(p, pgc) {}
    virtual ~sflow_amr_tsolver() = default;

    void start(lexer *p, ghostcell *pgc, slice &f, matrix2D &M, vec2D &xvec, vec2D &rhsvec, int var) override
    {
        auto around = [&](int a, int b)
        {
            const bool out = a < 0 || a >= p->knox || b < 0 || b >= p->knoy;
            return out && p->flagslice4[lij(p, a, b)] > 0;
        };

        if(var == 4 || var == 5)
        {
            int n = 0;
            SLICELOOP4
            {
                if(around(i + 1, j))
                {
                    rhsvec.V[n] -= M.n[n] * f(i + 1, j);
                    M.n[n] = 0.0;
                }
                if(around(i - 1, j))
                {
                    rhsvec.V[n] -= M.s[n] * f(i - 1, j);
                    M.s[n] = 0.0;
                }
                if(p->j_dir == 1)
                {
                    if(around(i, j + 1))
                    {
                        rhsvec.V[n] -= M.w[n] * f(i, j + 1);
                        M.w[n] = 0.0;
                    }
                    if(around(i, j - 1))
                    {
                        rhsvec.V[n] -= M.e[n] * f(i, j - 1);
                        M.e[n] = 0.0;
                    }
                }
                ++n;
            }
        }

        bs.start(p, pgc, f, M, xvec, rhsvec, var);
    }

private:
    sflow_bicgstab bs;
};

const int NTF = 3;  // k, epsilon/omega, eddy viscosity
}  // namespace

sflow_turb_io *sflow_amr::tio(int g)
{
    return static_cast<sflow_turb_io *>((g < 0) ? tb0 : SP(g)->pturb);
}

solver2D *sflow_amr::turb_solver(lexer *pp, ghostcell *pgc)
{
    return new sflow_amr_tsolver(pp, pgc);
}

void sflow_amr::turb_attach(sflow_turbulence *t)
{
    if(turb)
        tb0 = t;
}

void sflow_amr::turb_objects(sflow_amr_patch &c, ghostcell *pgc)
{
    if(!turb)
        return;

    if(p0->A260 == 1)
        c.pturb = new sflow_turb_ke_IM1(c.pp);
    else
        c.pturb = new sflow_turb_kw_IM1(c.pp);

    // the residual of the solver is per cell of the grid (the lexer of a patch has no global count)
    c.pp->cellnumtot2D = c.pp->knox * c.pp->knoy;

    c.ptsolv = turb_solver(c.pp, pgc);
}

// the face velocities and depths of a patch from its cell velocities, over all faces of its arrays (the production
// and the convection of the cells next to the box read the faces around them)
void sflow_amr::turb_faces(sflow_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm2D *b = c.b;

    for(int ii = pp->imin; ii < pp->imin + pp->imax - 1; ++ii)
        for(int jj = pp->jmin; jj < pp->jmin + pp->jmax; ++jj)
        {
            const bool w = pp->wet[lij(pp, ii, jj)] == 1 && pp->wet[lij(pp, ii + 1, jj)] == 1;
            b->P(ii, jj) = w ? 0.5 * (b->U(ii, jj) + b->U(ii + 1, jj)) : 0.0;
        }

    for(int ii = pp->imin; ii < pp->imin + pp->imax; ++ii)
        for(int jj = pp->jmin; jj < pp->jmin + pp->jmax - 1; ++jj)
        {
            const bool w = pp->wet[lij(pp, ii, jj)] == 1 && pp->wet[lij(pp, ii, jj + 1)] == 1;
            b->Q(ii, jj) = w ? 0.5 * (b->V(ii, jj) + b->V(ii, jj + 1)) : 0.0;
        }
}

// the cells around the level-l patches: from the patch of the level that holds them, else linear from the parent
void sflow_amr::turb_fill(ghostcell *pgc, int l)
{
    auto field = [&](int g, int n) -> slice &
    { return (n == 0) ? tio(g)->kin : (n == 1) ? tio(g)->eps
                                               : grid(g).b->eddyv; };

    fill_run(l, NTF, 8200 + l, [&](const reefamr_fill &f, double *v)
             {
                 for(int n=0; n<NTF; ++n)
                 {
                     if(f.kind==0)
                     v[n] = field(f.g,n)(f.si,f.sj);
                     else if(f.kind==1)
                     v[n] = sed_plin(field(f.g,n),grid(f.g).q,f.si,f.sj,f.ox,f.oy);
                     else
                     v[n] = 0.0;
                 } }, [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 for(int n=0; n<NTF; ++n)
                 field(id,n)(f.di,f.dj) = w[n]; });
}

void sflow_amr::turb_run(ghostcell *pgc, sflow_amr_patch &c)
{
    tscope ts(pgc, c.pp);
    c.pturb->start(c.pp, c.b, pgc, nullptr, c.pdiff, c.ptsolv, pflow_void);
}

void sflow_amr::turb_save(ghostcell *pgc, sflow_amr_patch &c)
{
    tscope ts(pgc, c.pp);
    c.pturb->ktimesave(c.pp, c.b, pgc);
    c.pturb->etimesave(c.pp, c.b, pgc);
}

// one turbulence step of the patches (sflow_f::mainloop, after the model of level 0)
void sflow_amr::turb_step(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(!turb || tb0 == nullptr || maxlev < 1 || patches_total == 0)
        return;

    {
        comms_off guard(pgc);
        for(auto q : P)
        {
            sflow_amr_patch *c = SP(q);
            c->pp->dt = p->dt;
            c->pp->dt_old = p->dt_old;
            c->pp->simtime = p->simtime;
            turb_faces(*c);
        }
    }

    // the patches, coarse to fine
    for(int l = 1; l <= maxlev; ++l)
    {
        turb_fill(pgc, l);
        for(int id : lev[l])
            turb_run(pgc, *SP(id));
    }

    // the covered cells: the mean of their children, finest first
    for(int l = maxlev; l >= 1; --l)
        block_up(l, NTF, 8260 + l, [&](reefamr_patch *q, int id, int k, double *v)
                 {
                 const int nby = q->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                 sflow_turb_io *t = tio(id);
                 fdm2D *bq = SP(q)->b;
                 for(int n=0; n<NTF; ++n)
                 {
                     slice &f = (n==0) ? t->kin : (n==1) ? t->eps : bq->eddyv;
                     v[n] = (sy==0) ? 0.5*(f(i0,j0) + f(i0+1,j0))
                                    : 0.25*(f(i0,j0) + f(i0+1,j0) + f(i0,j0+1) + f(i0+1,j0+1));
                 } }, [&](const reefamr_block &B, int key, const double *v)
                 {
                 if(std::isnan(v[0]))
                 return;
                 sflow_turb_io *t = tio(B.g);
                 fdm2D *bg = grid(B.g).b;
                 t->kin(B.ic,B.jc) = v[0];
                 t->eps(B.ic,B.jc) = v[1];
                 bg->eddyv(B.ic,B.jc) = v[2];
                 bg->kin(B.ic,B.jc) = v[0]; });

    // the ghost cells of level 0 (partition, walls; the labels of the models: k 20, epsilon/omega 30, eddy
    // viscosity 24)
    sflow_turb_io *t0 = tio(-1);
    pgc->gcsl_start4(p0, t0->kin, 20);
    pgc->gcsl_start4(p0, t0->eps, 30);
    pgc->gcsl_start4(p0, b0->kin, 20);
    pgc->gcsl_start4(p0, b0->eddyv, 24);

    // the old time level of the patches
    for(auto q : P)
        turb_save(pgc, *SP(q));
}

// regrid: the fields of the fresh level-l patches from the old patch of the level that held their cells, else linear
// from the parent (the blocks without an old patch, blkold of the state)
void sflow_amr::turb_regrid(ghostcell *pgc, int l, vector<reefamr_patch *> &oldP)
{
    if(!turb)
        return;

    {
        comms_off guard(pgc);
        for(int id : lev[l])
            if(P[id]->fresh)
            {
                sflow_amr_patch *c = SP(id);
                lexer *pp = c->pp;
                sflow_turb_io *t = tio(id);
                for(int ii = pp->imin; ii < pp->imin + pp->imax; ++ii)
                    for(int jj = pp->jmin; jj < pp->jmin + pp->jmax; ++jj)
                    {
                        t->kin(ii, jj) = 0.0;
                        t->eps(ii, jj) = 0.0;
                        c->b->eddyv(ii, jj) = 0.0;
                        c->b->kin(ii, jj) = 0.0;
                    }
            }
    }

    old_run(l, NTF, 8220 + l, oldP, [&](reefamr_patch *q, int io, int jo, double *v)
            {
                sflow_amr_patch *o = SP(q);
                sflow_turb_io *t = static_cast<sflow_turb_io*>(o->pturb);
                v[0] = t->kin(io,jo);
                v[1] = t->eps(io,jo);
                v[2] = o->b->eddyv(io,jo); }, [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
            {
                sflow_turb_io *t = tio(id);
                t->kin(ii,jj) = v[0];
                t->eps(ii,jj) = v[1];
                SP(q)->b->eddyv(ii,jj) = v[2];
                SP(q)->b->kin(ii,jj) = v[0]; });

    block_down_if(l, 4 * NTF, 8240 + l, [](reefamr_patch *q)
                  { return q->fresh; }, [&](const reefamr_block &B, int key, double *v)
                  {
                      gh G = grid(B.g);
                      sflow_turb_io *t = tio(B.g);
                      for(int n=0; n<NTF; ++n)
                      {
                          slice &f = (n==0) ? t->kin : (n==1) ? t->eps : G.b->eddyv;
                          for(int a=0; a<2; ++a)
                          for(int d=0; d<2; ++d)
                          v[4*n+2*a+d] = (d>sy) ? 0.0 : sed_plin(f,G.q,B.ic,B.jc,a==0?-1:1,(sy==0) ? 0 : (d==0?-1:1));
                      } }, [&](reefamr_patch *q, int id, int k, const double *v)
                  {
                      sflow_amr_patch *c = SP(q);
                      if(!q->fresh || c->blkold[k])
                      return;
                      const int nby = q->ny>>sy;
                      const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                      sflow_turb_io *t = tio(id);
                      for(int a=0; a<2; ++a)
                      for(int d=0; d<=sy; ++d)
                      {
                          t->kin(i0+a,j0+d) = v[2*a+d];
                          t->eps(i0+a,j0+d) = v[4+2*a+d];
                          c->b->eddyv(i0+a,j0+d) = v[8+2*a+d];
                          c->b->kin(i0+a,j0+d) = v[2*a+d];
                      } });

    // the old time level of the fresh patches
    for(int id : lev[l])
        if(P[id]->fresh)
            turb_save(pgc, *SP(id));
}
