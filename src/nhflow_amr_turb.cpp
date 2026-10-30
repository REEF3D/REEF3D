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
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_kepsilon_IM1.h"
#include"nhflow_komega_IM1.h"
#include"nhflow_scalar_ifou.h"
#include"nhflow_idiff.h"
#include"bicgstab_ijk.h"
#include"ioflow_amr_block.h"
#include"definitions.h"
#include<cmath>
#include<mpi.h>

using namespace nhflow_amr_detail;

// Turbulence (A 560 1/21 k-epsilon, 2/22 k-omega) on the grids of the mesh refinement and on multi-block grids.
//
// Every grid has its own turbulence model (nhflow_kepsilon_IM1 / nhflow_komega_IM1 on the lexer and the fdm of the
// grid) with its own scalar convection (first-order upwind), implicit diffusion and BiCGStab, as level 0. At the
// start of the step (before the momentum, as level 0 in driver_loop_nhflow):
//   1. the blocks 1..: the cells around the blocks take k, epsilon/omega and the eddy viscosity of the grid that holds
//      them (level 0 after its step, the other blocks from the end of the last step), every block steps its k and
//      epsilon/omega with them as boundary values; then the ring of the native grid from the blocks.  Four turned
//      blocks and one block agree to the solver tolerance (nhflow_multiblock/t2 with A 560: 5e-11 in k with N 43 and
//      N 44 1e-14); more passes do not help, a second step of a model from the same old time level is not the
//      same as the first (the sources take the new k and epsilon/omega)
//   2. levels 1..: the cells around the patches from a patch of the same level or linear from the parent (minmod
//      slopes: k and epsilon/omega stay positive; G 6: in the two halves of a layer), then the patches step
//   3. the covered cells take the mean of their children (per area on curvilinear grids), finest first; the cells
//      around the grids again (the eddy viscosity of both cells of a face is that of the grid that holds them)
//   4. the old time level of every grid (ktimesave, etimesave; level 0: the driver after the momentum)
// A fresh patch of a regrid takes its fields from an old patch of the same level where there was one, else linear
// from its parent.
//
// Subcycling (G 7 1, turb_level, turb_sync): level 0 makes its step in the driver as before; every step of a level
// l >= 1 (sub_level, with the time step of the level) starts with the turbulence step of its patches, the cells
// around them from a patch of the level or linear from the parent as it is then (the parent has made its own
// turbulence step at the start of its step: its fields are those of the end of that step, not interpolated in
// time, as level 0 for the patches without subcycling).  After the two steps of level l+1 the levels above l are
// restricted into the covered cells of level l and the cells around them filled again (sub_sync).  The time step
// of the model of a patch is that of its level.

namespace
{
typedef reefamr_comms_off comms_off;

inline double tmmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return fabs(a)<fabs(b) ? a : b;
}

const int NTF = 5;     // fields: k, epsilon/omega, EV, EV0, d->KIN
}

nhflow_amr::tfields nhflow_amr::tget(nhflow_turbulence *t)
{
    tfields f;
    if(nhflow_kepsilon_IM1 *a = dynamic_cast<nhflow_kepsilon_IM1*>(t))
    {
        f.k = a->KIN; f.e = a->EPS; f.kn = a->KN; f.en = a->EN;
    }
    else if(nhflow_komega_IM1 *b = dynamic_cast<nhflow_komega_IM1*>(t))
    {
        f.k = b->KIN; f.e = b->EPS; f.kn = b->KN; f.en = b->EN;
    }
    return f;
}

double* nhflow_amr::tfield(int g, int m)
{
    nhflow_turbulence *t = (g<0) ? tb0 : NP(g)->pturb;
    fdm_nhf *d = gfd(g);
    if(m==2) return d->EV;
    if(m==3) return d->EV0;
    if(m==4) return d->KIN;
    tfields f = tget(t);
    return (m==0) ? f.k : f.e;
}

void nhflow_amr::turb_attach(nhflow_turbulence *t)
{
    tb0 = t;
}

void nhflow_amr::turb_objects(nhflow_amr_patch &c, ghostcell *pgc)
{
    lexer *pp = c.pp;
    if(pp->A560==1 || pp->A560==21)
    c.pturb = new nhflow_kepsilon_IM1(pp,c.d,pgc);
    else
    c.pturb = new nhflow_komega_IM1(pp,c.d,pgc);
    c.ptconv = new nhflow_scalar_ifou(pp,1);
    c.ptdiff = new nhflow_idiff(pp);
    c.ptsolv = new bicgstab_ijk(pp,nullptr,pgc);
}

// the arrays of the model (the turbulence classes keep them for the whole run on level 0) and its helpers
void nhflow_amr::turb_free(nhflow_amr_patch &c)
{
    if(c.ptconv==nullptr)
    return;
    lexer *pp = c.pp;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);
    tfields f = tget(c.pturb);
    pp->del_Darray(f.k,n7);
    pp->del_Darray(f.e,n7);
    pp->del_Darray(f.kn,n7);
    pp->del_Darray(f.en,n7);
    if(nhflow_rans_io *r = dynamic_cast<nhflow_rans_io*>(c.pturb))
    {
        pp->del_Iarray(r->WALLF,n7);
        pp->del_Iarray(r->WETN,pp->imax*pp->jmax);
    }
    delete c.ptconv;
    delete c.ptdiff;
    delete static_cast<bicgstab_ijk*>(c.ptsolv);
    c.ptconv = nullptr; c.ptdiff = nullptr; c.ptsolv = nullptr;
}

void nhflow_amr::turb_run(ghostcell *pgc, nhflow_amr_patch &c)
{
    comms_off guard(pgc);
    pscope ps(pgc,c.d,d0);
    c.pturb->start(c.pp,c.d,pgc,c.ptconv,c.ptdiff,c.ptsolv,c.pflow,c.pvrans);
}

void nhflow_amr::turb_save(ghostcell *pgc, nhflow_amr_patch &c)
{
    comms_off guard(pgc);
    pscope ps(pgc,c.d,d0);
    c.pturb->ktimesave(c.pp,c.d,pgc);
    c.pturb->etimesave(c.pp,c.d,pgc);
}

// the wall and partition ghost cells of the fields of grid g (as at the end of the model's step)
void nhflow_amr::turb_bc(ghostcell *pgc, int g)
{
    lexer *q = glex(g);
    fdm_nhf *d = gfd(g);
    double *k = tfield(g,0), *e = tfield(g,1);
    if(g<0)
    {
        pgc->start20V(q,k,20);
        pgc->start30V(q,e,30);
        pgc->start24V(q,d->EV,24);
        pgc->start24V(q,d->EV0,24);
        pgc->start4V(q,d->KIN,20);
        return;
    }
    comms_off guard(pgc);
    pscope ps(pgc,d,d0);
    pgc->start20V(q,k,20);
    pgc->start30V(q,e,30);
    pgc->start24V(q,d->EV,24);
    pgc->start24V(q,d->EV0,24);
    pgc->start4V(q,d->KIN,20);
}

// level 0: the cells around the blocks and the ring of the native grid from the grid that holds them
void nhflow_amr::turb_fill0(ghostcell *pgc)
{
    const int K = p0->knoz;
    const int nv = NTF*K;

    auto put = [&](int g, const reefamr_fill &f, const double *w)
    {
        lexer *q = glex(g);
        for(int m=0; m<NTF; ++m)
        {
            double *a = tfield(g,m);
            for(int k=0; k<K; ++k)
            a[cidx(q,f.di,f.dj,k)] = w[m*K+k];
        }
    };

    fill_run0(nv,7095,
              [&](const reefamr_fill &f, double *v)
              {
                  if(f.kind!=0)
                  {
                      for(int m=0; m<nv; ++m)
                      v[m] = 0.0;
                      return;
                  }
                  lexer *q = glex(f.g);
                  for(int m=0; m<NTF; ++m)
                  {
                      const double *a = tfield(f.g,m);
                      for(int k=0; k<K; ++k)
                      v[m*K+k] = a[cidx(q,f.si,f.sj,k)];
                  }
              },
              [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w) { put(id,f,w); },
              [&](const reefamr_fill &f, const double *w) { put(-1,f,w); });

    // the ring cells in the partition halos of the native grid
    if(p0->mpi_size>1)
    {
        for(int m=0; m<NTF; ++m)
        {
            double *a = tfield(-1,m);
            pgc->gcparaxV(p0,a,20);
            pgc->gcparacoxV(p0,a,20);
        }
    }
}

// levels 1..: the cells around the level-l patches, from a patch of the same level or linear from the parent
void nhflow_amr::turb_fill(ghostcell *pgc, int l)
{
    const int K = klev(l);
    const int Kc = klev(l-1);
    const int nv = NTF*K;
    vector<double> uc(Kc);

    fill_run(l,nv,7096+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 {
                     lexer *q = glex(f.g);
                     for(int m=0; m<NTF; ++m)
                     {
                         const double *a = tfield(f.g,m);
                         for(int k=0; k<K; ++k)
                         v[m*K+k] = a[cidx(q,f.si,f.sj,k)];
                     }
                     return;
                 }
                 if(f.kind==1)
                 {
                     lexer *q = glex(f.g);
                     auto ok = [&](int a, int b) { const int n = lij(q,a,b); return q->flagslice4[n]>0 && q->wet[n]==1; };
                     const bool ex = ok(f.si+1,f.sj) && ok(f.si-1,f.sj);
                     const bool ey = ok(f.si,f.sj+1) && ok(f.si,f.sj-1);
                     for(int m=0; m<NTF; ++m)
                     {
                         const double *a = tfield(f.g,m);
                         for(int k=0; k<Kc; ++k)
                         {
                             const double c = a[cidx(q,f.si,f.sj,k)];
                             double sx = 0.0, sy = 0.0;
                             if(ex)
                             sx = tmmod(a[cidx(q,f.si+1,f.sj,k)]-c,c-a[cidx(q,f.si-1,f.sj,k)]);
                             if(ey)
                             sy = tmmod(a[cidx(q,f.si,f.sj+1,k)]-c,c-a[cidx(q,f.si,f.sj-1,k)]);
                             uc[k] = c + 0.25*f.ox*sx + 0.25*f.oy*sy;
                         }
                         vcell(q,&uc[0],Kc,K/Kc,&v[m*K]);
                     }
                     return;
                 }
                 for(int m=0; m<nv; ++m)
                 v[m] = 0.0;
             },
             [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
             {
                 lexer *pp = q->pp;
                 for(int m=0; m<NTF; ++m)
                 {
                     double *a = tfield(id,m);
                     for(int k=0; k<K; ++k)
                     a[cidx(pp,f.di,f.dj,k)] = (m<2 || m==4) ? MAX(w[m*K+k],0.0) : w[m*K+k];
                 }
             });
}

// the covered cells: the mean of their children (per area on curvilinear grids; G 6: of both fine layers)
void nhflow_amr::turb_restrict()
{
    for(int l=wtop(); l>=wlo+1; --l)
    {
        const int K = klev(l-1);
        block_up(l,NTF*K,7340+l,
                 [&](reefamr_patch *q, int id, int k, double *v)
                 {
                     lexer *pp = q->pp;
                     const int nby = q->ny>>sy;
                     const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                     const int fz = pp->knoz/K;
                     double w[4];
                     cw4(pp,i0,j0,w);
                     for(int m=0; m<NTF; ++m)
                     {
                         const double *f = tfield(id,m);
                         for(int kk=0; kk<K; ++kk)
                         {
                             double r = 0.0;
                             for(int kf=fz*kk; kf<fz*kk+fz; ++kf)
                             r += (sy==0 ? w[0]*f[cidx(pp,i0,j0,kf)]+w[1]*f[cidx(pp,i0+1,j0,kf)] : w[0]*f[cidx(pp,i0,j0,kf)]+w[1]*f[cidx(pp,i0+1,j0,kf)]+w[2]*f[cidx(pp,i0,j0+1,kf)]+w[3]*f[cidx(pp,i0+1,j0+1,kf)]);
                             v[m*K+kk] = r/double(fz);
                         }
                     }
                 },
                 [&](const reefamr_block &B, int key, const double *v)
                 {
                     lexer *q = glex(B.g);
                     for(int m=0; m<NTF; ++m)
                     {
                         double *f = tfield(B.g,m);
                         for(int kk=0; kk<K; ++kk)
                         f[cidx(q,B.ic,B.jc,kk)] = v[m*K+kk];
                     }
                 });
    }
}

// the blocks 1.. start as level 0 (nhflow_rans_io::ini: the inflow profile with B 60), then the old time level
void nhflow_amr::turb_blocks_ini(ghostcell *pgc)
{
    if(!turb)
    return;
    for(int id : lev[0])
    {
        nhflow_amr_patch *c = NP(id);
        lexer *pp = c->pp;
        comms_off guard(pgc);
        pscope ps(pgc,c->d,d0);
        // the inflow profile in all cells as on level 0 (a block without an inflow face has B 60 0), with the
        // inflow velocity of level 0 where the block has no inflow of its own
        const int b60 = pp->B60;
        const double ui = pp->Ui;
        if(p0->B60>=1)
        {
            pp->B60 = p0->B60;
            if(c->pio==nullptr || c->pio->in==0)
            pp->Ui = p0->Ui;
        }
        c->pturb->ini(pp,c->d,pgc);
        pp->B60 = b60;
        pp->Ui = ui;
    }
    turb_fill0(pgc);
    for(int id : lev[0])
    turb_save(pgc,*NP(id));
}

// the fresh patches of level l after a regrid: linear from the parent, then from an old patch of the same level
void nhflow_amr::turb_regrid(ghostcell *pgc, int l, vector<reefamr_patch*> &oldP)
{
    if(!turb)
    return;

    const int K = klev(l);
    const int Kc = klev(l-1);

    for(int id : lev[l])
    if(P[id]->fresh)
    {
        nhflow_amr_patch *c = NP(id);
        comms_off guard(pgc);
        pscope ps(pgc,c->d,d0);
        c->pturb->ini(c->pp,c->d,pgc);
    }

    // the 2x2 children of every parent of a fresh patch: per child and field the K layers
    block_down_if(l,4*NTF*K,7640+l,[](reefamr_patch *q) { return q->fresh; },
                  [&](const reefamr_block &B, int key, double *v)
                  {
                      lexer *q = glex(B.g);
                      auto ok = [&](int a, int b) { const int n = lij(q,a,b); return q->flagslice4[n]>0 && q->wet[n]==1; };
                      const bool ex = ok(B.ic+1,B.jc) && ok(B.ic-1,B.jc);
                      const bool ey = ok(B.ic,B.jc+1) && ok(B.ic,B.jc-1);
                      vector<double> uc(Kc);
                      for(int a=0; a<2; ++a)
                      for(int b=0; b<=sy; ++b)
                      for(int m=0; m<NTF; ++m)
                      {
                          const double *f = tfield(B.g,m);
                          for(int k=0; k<Kc; ++k)
                          {
                              const double c = f[cidx(q,B.ic,B.jc,k)];
                              double sx = 0.0, sy = 0.0;
                              if(ex)
                              sx = tmmod(f[cidx(q,B.ic+1,B.jc,k)]-c,c-f[cidx(q,B.ic-1,B.jc,k)]);
                              if(ey)
                              sy = tmmod(f[cidx(q,B.ic,B.jc+1,k)]-c,c-f[cidx(q,B.ic,B.jc-1,k)]);
                              uc[k] = c + 0.25*(a==0?-1:1)*sx + 0.25*(b==0?-1:1)*sy;
                          }
                          vcell(q,&uc[0],Kc,K/Kc,&v[((2*a+b)*NTF+m)*K]);
                      }
                  },
                  [&](reefamr_patch *q, int id, int k, const double *v)
                  {
                      if(!q->fresh)
                      return;
                      lexer *pp = q->pp;
                      const int nby = q->ny>>sy;
                      const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                      for(int a=0; a<2; ++a)
                      for(int b=0; b<=sy; ++b)
                      for(int m=0; m<NTF; ++m)
                      {
                          double *f = tfield(id,m);
                          for(int kk=0; kk<K; ++kk)
                          f[cidx(pp,i0+a,j0+b,kk)] = v[((2*a+b)*NTF+m)*K+kk];
                      }
                  });

    // where an old patch of the same level was, its fields
    old_run(l,NTF*K,7650+l,oldP,
            [&](reefamr_patch *o, int io, int jo, double *v)
            {
                nhflow_amr_patch *oc = NP(o);
                lexer *op = o->pp;
                tfields tf = tget(oc->pturb);
                const double *src[NTF] = {tf.k,tf.e,oc->d->EV,oc->d->EV0,oc->d->KIN};
                for(int m=0; m<NTF; ++m)
                for(int kk=0; kk<K; ++kk)
                v[m*K+kk] = src[m][cidx(op,io,jo,kk)];
            },
            [&](reefamr_patch *c, int id, int ii, int jj, const double *v)
            {
                lexer *pp = c->pp;
                for(int m=0; m<NTF; ++m)
                {
                    double *f = tfield(id,m);
                    for(int kk=0; kk<K; ++kk)
                    f[cidx(pp,ii,jj,kk)] = v[m*K+kk];
                }
            });

    turb_fill(pgc,l);

    for(int id : lev[l])
    if(P[id]->fresh)
    {
        turb_bc(pgc,id);
        turb_save(pgc,*NP(id));
    }
}

// one turbulence step of all grids (nhflow_amr::step, before the momentum; level 0 has made its own)
void nhflow_amr::turb_step(lexer *p, ghostcell *pgc)
{
    if(!turb)
    return;

    const double t0 = MPI_Wtime();

    // the blocks (level 0 has made its step: the cells around the blocks from it and from the blocks at the end of
    // the last step), then the ring of the native grid from the blocks
    if(nblk>0)
    {
        turb_fill0(pgc);
        for(int id : lev[0])
        turb_run(pgc,*NP(id));
        turb_fill0(pgc);
    }

    // the patches, coarse to fine
    for(int l=1; l<=maxlev; ++l)
    {
        turb_fill(pgc,l);
        for(int id : lev[l])
        turb_run(pgc,*NP(id));
    }

    // the covered cells, then the cells around the grids and the ghost cells of level 0
    if(maxlev>=1)
    {
        {
        comms_off guard(pgc);
        turb_restrict();
        }
        turb_bc(pgc,-1);
        if(nblk>0)
        turb_fill0(pgc);
        for(int l=1; l<=maxlev; ++l)
        turb_fill(pgc,l);
    }

    for(auto q : P)
    turb_save(pgc,*NP(q));

    tturb += MPI_Wtime()-t0;
    ++tturb_n;
}

// G 7 1: the turbulence step of the level-l patches at the start of their step (sub_level): the cells around them,
// the step with the time step of the level, the old time level
void nhflow_amr::turb_level(ghostcell *pgc, int l)
{
    if(!turb)
    return;

    const double t0 = MPI_Wtime();

    turb_fill(pgc,l);
    for(int id : lev[l])
    turb_run(pgc,*NP(id));
    for(int id : lev[l])
    turb_save(pgc,*NP(id));

    tturb += MPI_Wtime()-t0;
    ++tturb_n;
}

// G 7 1: after the steps of the levels above l (sub_sync): their fields restricted into the covered cells down to
// level l, the ghost cells of level 0, the cells around the patches above l again; the old time level of the level-l
// patches with their restricted cells (level 0: the driver after the momentum)
void nhflow_amr::turb_sync(ghostcell *pgc, int l)
{
    if(!turb || l>=maxlev)
    return;

    const double t0 = MPI_Wtime();

    wlo = l;
    whi = maxlev;
    {
    comms_off guard(pgc);
    turb_restrict();
    }
    wlo = 0;
    whi = -1;

    if(l==0)
    turb_bc(pgc,-1);
    else
    for(int id : lev[l])
    turb_bc(pgc,id);

    for(int ll=MAX(l,1); ll<=maxlev; ++ll)
    turb_fill(pgc,ll);

    if(l>=1)
    for(int id : lev[l])
    turb_save(pgc,*NP(id));

    tturb += MPI_Wtime()-t0;
}
