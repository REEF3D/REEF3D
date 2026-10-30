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

#include"reefamr.h"
#include"reef3d_mpmd.h"
#include"reefamr_krylov.h"
#include"reefamr_xzlu.h"
#include"reefmg_core.h"
#include"reefmg2D.h"
#include"si_eta2D.h"
#include"matrix2D.h"
#include"vec2D.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice4.h"
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<map>
#include<cmath>
#include<limits>
#include<mpi.h>

//  The semi-implicit free surface on the grids of the mesh refinement: one eta system over the leaf cells of all
//  grids (see reefamr.h).  The level-0 part (the native grid and the blocks 1.. of a multi-block grid: their own
//  rows, the cells around them filled from the cells they stand for, restricted additive Schwarz with one REEFMG2D
//  V-cycle per grid, G 57 coarse space) is the solve of nhflow_amr_si.cpp, moved here unchanged; the levels add the
//  covered cells (no row, the mean of their children), the interpolated cells around the patches, the fine fluxes on
//  the coarse faces next to the patches and the FAC sweep of the preconditioner (as sflow_amr_nh.cpp).

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

typedef reefamr_comms_off comms_off;

// the vectors of a grid
enum { VB=0, VR=1, VRH=2, VPV=3, VVV=4, VS=5, VT=6, VPH=7, VSH=8, VX=9, VZC=10, VAZ=11, VRT=12, VAG=13, VNV=14 };

// passes of the patch-local corrections per level (as sflow_amr_nh.cpp)
const int RSPASS = 2;
}

// one grid of the eta system
struct reefamr::rsi_g
{
    int g = -1, lev = 0;
    lexer *q = nullptr;
    rsi_rows R;
    vector<int> own;            // own cells: i, j, row of the si_eta2D assembly (covered cells included)
    vector<int> leaf;           // the leaf cells among them
    vector<int> kind;           // per cell (lij): -2 no row, -1 covered, row n of a leaf cell
    vector<int> rowof;          // per cell (lij): row of an own cell, -1 none
    int nrow = 0;
    vector<slice*> v;
    // level 0: the V-cycle of the grid (restricted additive Schwarz)
    reefmg2D *mg = nullptr;
    vec2D *rv = nullptr;
    slice *z = nullptr;
    char rebuild = 1;
    // levels: the patch multigrid (MPI_COMM_SELF), on a 2D grid the exact line solve
    reefmg_core *pmg = nullptr;
    reefamr_xzlu xz;
    // the a-fluxes a_f (x_R - x_L) of the box faces of a patch (operator), per side
    vector<double> fx[4];

    ~rsi_g()
    {
        for(auto s : v)
        delete s;
        delete mg;
        delete rv;
        delete z;
        delete pmg;
    }
};

struct reefamr::rsi_csi
{
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    vector<double> rl, rc;
};

namespace
{
// the vector space of the composite eta solve for reefamr_bicgstab
struct rsi_space
{
    reefamr *a;
    static const int B=VB, R=VR, RH=VRH, PV=VPV, VV=VVV, S=VS, T=VT, PH=VPH, SH=VSH;
    void apply(int x, int y) { a->rsi_apply(x,y); }
    void prec(int r, int z) { a->rsi_prec(r,z); }
    double dot(int x, int y) { return a->rsi_dot(x,y); }
    void op_start() { a->rsi_opstart(); }
    void op_p(double beta, double om) { a->rsi_op_p(beta,om); }
    void op_s(double alp) { a->rsi_op_s(alp); }
    void op_x(double alp, double om) { a->rsi_op_x(alp,om); }
    void op_restart() { a->rsi_restart(); }
};
}

// ------------------------------------------------------------------------------------------- grids
void reefamr::rsi_run(int g, const std::function<void()> &f)
{
    if(g<0)
    {
        f();
        return;
    }

    comms_off guard(pgc0);
    f();
}

bool reefamr::rsi_own(int g, int ii, int jj)
{
    lexer *q = glex(g);

    if(q->flagslice4[lij(q,ii,jj)]<0)
    return false;

    if(g<0)
    return !(nblk>0 && nring[lij(q,ii,jj)]==1);

    reefamr_patch *c = P[g];
    return ii>=EXT && ii<EXT+c->nx && jj>=EY && jj<EY+c->ny;
}

bool reefamr::rsi_leaf(int g, int ii, int jj)
{
    if(!rsi_own(g,ii,jj))
    return false;

    const int l = (g<0) ? 0 : P[g]->lev;
    if(l>=maxlev)
    return true;

    const int I = (g<0) ? ii+O0i : ii-EXT+P[g]->I0;
    const int J = (g<0) ? jj+O0j : jj-EY+P[g]->J0;

    return !covered(l+1,2*I,ychild(J));
}

const vector<int>& reefamr::rsi_grids()
{
    if(rsi_layout!=regrids || rsg.empty())
    rsi_setup();

    return rsig;
}

void reefamr::rsi_free()
{
    for(auto c : rsg)
    delete c;
    rsg.clear();
    rsig.clear();
}

void reefamr::rsi_setup()
{
    rsi_free();
    rsi_layout = regrids;

    rsig.push_back(-1);
    for(int id : lev[0])
    rsig.push_back(id);

    rsi_lev = false;
    for(int l=1; l<=maxlev && l<(int)lev.size(); ++l)
    for(int id : lev[l])
    {
        rsig.push_back(id);
        rsi_lev = true;
    }

    for(int g : rsig)
    {
        rsi_g *c = new rsi_g;
        c->g = g;
        c->lev = (g<0) ? 0 : P[g]->lev;
        c->q = glex(g);
        rsi_rows_of(g,c->R);
        lexer *q = c->q;

        c->kind.assign((size_t)q->imax*q->jmax,-2);
        c->rowof.assign((size_t)q->imax*q->jmax,-1);

        // the rows in the order of SLICELOOP4 (si_eta2D::build)
        int n=0;
        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        if(q->flagslice4[lij(q,ii,jj)]>0)
        {
            if(rsi_own(g,ii,jj))
            {
                c->own.push_back(ii);
                c->own.push_back(jj);
                c->own.push_back(n);
                c->rowof[lij(q,ii,jj)] = n;

                if(rsi_leaf(g,ii,jj))
                {
                    c->leaf.push_back(ii);
                    c->leaf.push_back(jj);
                    c->leaf.push_back(n);
                    c->kind[lij(q,ii,jj)] = n;
                }
                else
                c->kind[lij(q,ii,jj)] = -1;
            }
        ++n;
        }

        c->nrow = n;

        for(int v=0; v<VNV; ++v)
        {
            slice *s = new slice4(q);
            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
            (*s)(ii,jj) = 0.0;
            c->v.push_back(s);
        }

        // level 0: the V-cycle preconditioner of every grid (native: over its ranks; blocks: on their rank)
        if(c->lev==0)
        {
            if(g<0)
            c->mg = new reefmg2D(q,rsipgc);

            else
            rsi_run(g,[&](){ c->mg = new reefmg2D(q,rsipgc,true); });

            c->rv = new vec2D(q,max(c->nrow,1));
            c->z = new slice4(q);
        }

        if(c->lev>0)
        {
            reefamr_patch *P0 = P[g];
            c->fx[0].assign(P0->ny,0.0);
            c->fx[1].assign(P0->ny,0.0);
            c->fx[2].assign(sy==1 ? P0->nx : 0,0.0);
            c->fx[3].assign(sy==1 ? P0->nx : 0,0.0);
        }

        rsg.push_back(c);
    }

    // the coarse space (G 57) and the rows of the patch multigrids are made per solve
}

// ------------------------------------------------------------------------------------------- fills
// covered cells: the mean of their children (curvilinear grid: area weighted), finest first; the parent may lie on
// another rank (block plans)
void reefamr::rsi_restrict(int tag, const std::function<slice&(int)> &F)
{
    for(int l=maxlev; l>=1; --l)
    block_up(l,1,tag+l,
             [&](reefamr_patch *c, int id, int k, double *v)
             {
                 lexer *pp = c->pp;
                 slice &x = F(id);
                 const int nby = c->ny>>sy;
                 const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);

                 if(pp->geo_curv==1)
                 {
                     double s=0.0, a=0.0;
                     for(int u=0; u<2; ++u)
                     for(int d=0; d<=sy; ++d)
                     {
                         const double A = pp->AREA2D[lij(pp,i0+u,j0+d)];
                         s += A*x(i0+u,j0+d);
                         a += A;
                     }
                     v[0] = a>0.0 ? s/a : 0.0;
                     return;
                 }

                 v[0] = (sy==0) ? 0.5*(x(i0,j0)+x(i0+1,j0)) : 0.25*(x(i0,j0)+x(i0+1,j0)+x(i0,j0+1)+x(i0+1,j0+1));
             },
             [&](const reefamr_block &B, int key, const double *v)
             {
                 F(B.g)(B.ic,B.jc) = v[0];
             });
}

// linear interpolation of the coarse cell (ic,jc) of a grid into its child in quadrant (ox,oy) (central slopes,
// none across a wall: as sflow_amr::nh_eval)
double reefamr::rsi_interp(lexer *q, slice &x, int ic, int jc, int ox, int oy)
{
    const double x0 = x(ic,jc);

    auto val = [&](int a, int b) { return q->flagslice4[lij(q,a,b)]<0 ? x0 : x(a,b); };

    return x0 + 0.125*(val(ic+1,jc)-val(ic-1,jc))*ox + 0.125*(val(ic,jc+1)-val(ic,jc-1))*oy;
}

void reefamr::rsi_fill(int tag, const std::function<slice&(int)> &F)
{
    // the cells around the level-0 grids, the ring and the halo of the native grid (without patches the fill of the
    // multi-block solve before this layer)
    if(!rsi_lev)
    {
        xslice(1,tag,[&](int g, int k) -> slice* { return &F(g); });
        return;
    }

    rsi_restrict(tag+10,F);

    const double nan = std::numeric_limits<double>::quiet_NaN();

    if(nblk>0)
    fill_run0(1,tag,
              [&](const reefamr_fill &f, double *v) { v[0] = (f.kind==0) ? F(f.g)(f.si,f.sj) : nan; },
              [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w) { if(w[0]==w[0]) F(id)(f.di,f.dj) = w[0]; },
              [&](const reefamr_fill &f, const double *w) { if(w[0]==w[0]) F(-1)(f.di,f.dj) = w[0]; });

    native_gc4(F(-1));

    // the cells around the patches, coarse to fine: a patch of the same level, or the coarser level interpolated
    for(int l=1; l<=maxlev; ++l)
    fill_run(l,1,tag+20+l,
             [&](const reefamr_fill &f, double *v)
             {
                 if(f.kind==0)
                 v[0] = F(f.g)(f.si,f.sj);
                 else if(f.kind==1)
                 v[0] = rsi_interp(glex(f.g),F(f.g),f.si,f.sj,f.ox,f.oy);
                 else
                 v[0] = 0.0;
             },
             [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w) { F(id)(f.di,f.dj) = w[0]; });
}

// Krylov vector v of all grids
void reefamr::rsi_fillv(int v, bool all)
{
    auto F = [&](int g) -> slice&
    {
        for(auto c : rsg)
        if(c->g==g)
        return *c->v[v];
        return *rsg[0]->v[v];
    };

    if(all)
    {
        rsi_fill(7620,F);
        return;
    }

    // the level-0 grids only (the level-0 part of the preconditioner)
    xslice(1,7620,[&](int g, int k) -> slice* { return &F(g); });
}

// ------------------------------------------------------------------------------------------- operator
// y = A x on the own rows of the level-0 grids, the cells around them from the level-0 fill (no patch values)
void reefamr::rsi_rows0(int x, int y)
{
    rsi_fillv(x,!rsi_lev ? true : false);

    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;

        const matrix2D &M = c->R.e->mat();
        slice &X = *c->v[x];
        slice &Y = *c->v[y];
        const vector<int> &r = c->own;

        for(size_t q=0; q<r.size(); q+=3)
        {
            const int ii=r[q], jj=r[q+1], n=r[q+2];

            Y(ii,jj) = M.p[n]*X(ii,jj) + M.n[n]*X(ii+1,jj) + M.s[n]*X(ii-1,jj)
                     + M.w[n]*X(ii,jj+1) + M.e[n]*X(ii,jj-1);
        }
    }
}

void reefamr::rsi_apply(int x, int y)
{
    if(!rsi_lev)
    {
        rsi_rows0(x,y);
        return;
    }

    rsi_fillv(x,true);

    // the leaf rows; the covered rows are 0
    for(auto c : rsg)
    {
        const matrix2D &M = c->R.e->mat();
        slice &X = *c->v[x];
        slice &Y = *c->v[y];
        const vector<int> &r = c->leaf;

        for(size_t q=0; q<c->own.size(); q+=3)
        Y(c->own[q],c->own[q+1]) = 0.0;

        for(size_t q=0; q<r.size(); q+=3)
        {
            const int ii=r[q], jj=r[q+1], n=r[q+2];

            Y(ii,jj) = M.p[n]*X(ii,jj) + M.n[n]*X(ii+1,jj) + M.s[n]*X(ii-1,jj)
                     + M.w[n]*X(ii,jj+1) + M.e[n]*X(ii,jj-1);
        }
    }

    // the coarse faces next to the patches: the mean of the fine a-fluxes a_f (x_R - x_L) instead of the coarse one
    rsm_begin(1);

    for(auto c : rsg)
    {
        if(c->lev==0)
        continue;

        reefamr_patch *P0 = P[c->g];
        slice &X = *c->v[x];
        slice &ax = *c->R.ax;
        slice &ay = *c->R.ay;
        const int il = EXT-1, ih = EXT+P0->nx-1, jl = EY-1, jh = EY+P0->ny-1;

        for(int r=0; r<P0->ny; ++r)
        {
            const int jj = EY+r;
            *rsm_slot(c->g,0,r) = ax(il,jj)*(X(il+1,jj) - X(il,jj));
            *rsm_slot(c->g,1,r) = ax(ih,jj)*(X(ih+1,jj) - X(ih,jj));
        }

        if(sy==1)
        for(int r=0; r<P0->nx; ++r)
        {
            const int ii = EXT+r;
            *rsm_slot(c->g,2,r) = ay(ii,jl)*(X(ii,jl+1) - X(ii,jl));
            *rsm_slot(c->g,3,r) = ay(ii,jh)*(X(ii,jh+1) - X(ii,jh));
        }
    }

    for(int l=maxlev; l>=1; --l)
    rsm_send(l,8600);

    for(auto c : rsg)
    {
        lexer *q = c->q;
        const matrix2D &M = c->R.e->mat();
        slice &X = *c->v[x];
        slice &Y = *c->v[y];
        const bool cv = (q->geo_curv==1);

        // the geometry factor of a face of a cell (the face length over the cell area)
        auto geo = [&](int dir, int fi, int fj, int ci, int cj)
        {
            if(cv)
            return (dir==1 ? q->LEN1[lij(q,fi,fj)] : q->LEN2[lij(q,fi,fj)])/q->AREA2D[lij(q,ci,cj)];
            return dir==1 ? 1.0/q->DXN[ci+marge] : 1.0/q->DYN[cj+marge];
        };

        auto isolated = [&](int n) { return M.n[n]==0.0 && M.s[n]==0.0 && M.w[n]==0.0 && M.e[n]==0.0; };

        rsm_take(c->g,[&](int dir, int fi, int fj, const double *v, int sgn)
        {
            const double F = sgn*v[0];
            const int ai = fi, aj = fj;                         // the cell on the low side of the face
            const int bi = dir==1 ? fi+1 : fi, bj = dir==1 ? fj : fj+1;   // ... on the high side
            const int ka = c->kind[lij(q,ai,aj)];
            const int kb = c->kind[lij(q,bi,bj)];

            if(ka>=0 && !isolated(ka))
            {
                const double co = (dir==1) ? M.n[ka] : M.w[ka];
                Y(ai,aj) += -geo(dir,fi,fj,ai,aj)*F - co*(X(bi,bj) - X(ai,aj));
            }

            if(kb>=0 && !isolated(kb))
            {
                const double co = (dir==1) ? M.s[kb] : M.e[kb];
                Y(bi,bj) += geo(dir,fi,fj,bi,bj)*F - co*(X(ai,aj) - X(bi,bj));
            }
        });
    }
}

double reefamr::rsi_dot(int a, int b)
{
    double s=0.0;

    for(auto c : rsg)
    {
        slice &A = *c->v[a];
        slice &B = *c->v[b];
        const vector<int> &w = c->leaf;

        for(size_t q=0; q<w.size(); q+=3)
        s += A(w[q],w[q+1])*B(w[q],w[q+1]);
    }

    double r=0.0;
    MPI_Allreduce(&s,&r,1,MPI_DOUBLE,MPI_SUM,rsipgc->mpi_comm);
    return r;
}

// the vector operations of reefamr_bicgstab over the leaf cells
void reefamr::rsi_opstart()
{
    for(auto c : rsg)
    {
        const vector<int> &w = c->leaf;
        for(size_t q=0; q<w.size(); q+=3)
        {
            const int ii=w[q], jj=w[q+1];
            const double r = (*c->v[VB])(ii,jj) - (*c->v[VVV])(ii,jj);
            (*c->v[VR])(ii,jj) = r;
            (*c->v[VRH])(ii,jj) = r;
            (*c->v[VPV])(ii,jj) = 0.0;
            (*c->v[VVV])(ii,jj) = 0.0;
        }
    }
}

void reefamr::rsi_op_p(double beta, double om)
{
    for(auto c : rsg)
    {
        const vector<int> &w = c->leaf;
        for(size_t q=0; q<w.size(); q+=3)
        {
            const int ii=w[q], jj=w[q+1];
            (*c->v[VPV])(ii,jj) = (*c->v[VR])(ii,jj) + beta*((*c->v[VPV])(ii,jj) - om*(*c->v[VVV])(ii,jj));
        }
    }
}

void reefamr::rsi_op_s(double alp)
{
    for(auto c : rsg)
    {
        const vector<int> &w = c->leaf;
        for(size_t q=0; q<w.size(); q+=3)
        {
            const int ii=w[q], jj=w[q+1];
            (*c->v[VS])(ii,jj) = (*c->v[VR])(ii,jj) - alp*(*c->v[VVV])(ii,jj);
        }
    }
}

void reefamr::rsi_op_x(double alp, double om)
{
    for(auto c : rsg)
    {
        const vector<int> &w = c->leaf;
        for(size_t q=0; q<w.size(); q+=3)
        {
            const int ii=w[q], jj=w[q+1];
            (*c->v[VX])(ii,jj) += alp*(*c->v[VPH])(ii,jj) + om*(*c->v[VSH])(ii,jj);
            (*c->v[VR])(ii,jj) = (*c->v[VS])(ii,jj) - om*(*c->v[VT])(ii,jj);
        }
    }
}

void reefamr::rsi_restart()
{
    for(auto c : rsg)
    {
        const vector<int> &w = c->leaf;
        for(size_t q=0; q<w.size(); q+=3)
        {
            const int ii=w[q], jj=w[q+1];
            (*c->v[VRH])(ii,jj) = (*c->v[VR])(ii,jj);
            (*c->v[VPV])(ii,jj) = 0.0;
            (*c->v[VVV])(ii,jj) = 0.0;
        }
    }
}

// ------------------------------------------------------------------------------------------- preconditioner
// G 57: the two-level preconditioner of the level-0 grids, the coarse correction zc = P Ac^-1 R r first, then the
// V-cycles of the grids on the remaining residual r - A zc: z = zc + RAS(r - A zc); without G 57 the V-cycles alone;
// with patches the FAC sweep over the levels follows
void reefamr::rsi_prec(int r, int z)
{
    if(rsi_lev)
    rsi_restrict(8610,[&](int g) -> slice& { for(auto c : rsg) if(c->g==g) return *c->v[r]; return *rsg[0]->v[r]; });

    if(rsics_n<=0)
    rsi_ras(r,z);

    else
    {
        rsi_cs_correct(r,VZC);
        rsi_rows0(VZC,VAZ);

        for(auto c : rsg)
        {
            if(c->lev>0)
            continue;
            slice &R = *c->v[r];
            slice &RT = *c->v[VRT];
            slice &AZ = *c->v[VAZ];
            const vector<int> &w = c->own;
            for(size_t k=0; k<w.size(); k+=3)
            RT(w[k],w[k+1]) = R(w[k],w[k+1]) - AZ(w[k],w[k+1]);
        }

        rsi_ras(VRT,z);

        for(auto c : rsg)
        {
            if(c->lev>0)
            continue;
            slice &Z = *c->v[z];
            slice &ZC = *c->v[VZC];
            const vector<int> &w = c->own;
            for(size_t k=0; k<w.size(); k+=3)
            Z(w[k],w[k+1]) += ZC(w[k],w[k+1]);
        }
    }

    if(!rsi_lev)
    return;

    // FAC: level by level the correction interpolated into the patches, then the patch-local corrections of the
    // residual it leaves (covered cells included: they carry the correction of the next finer level)
    auto Zf = [&](int g) -> slice& { for(auto c : rsg) if(c->g==g) return *c->v[z]; return *rsg[0]->v[z]; };

    // the level-0 correction around the level-0 grids
    xslice(1,8620,[&](int g, int k) -> slice* { return &Zf(g); });

    for(int l=1; l<=maxlev; ++l)
    {
        rsi_prolong(l,z);

        auto fill_l = [&]()
        {
            fill_run(l,1,8640+l,
                     [&](const reefamr_fill &f, double *v)
                     {
                         if(f.kind==0)
                         v[0] = Zf(f.g)(f.si,f.sj);
                         else if(f.kind==1)
                         v[0] = rsi_interp(glex(f.g),Zf(f.g),f.si,f.sj,f.ox,f.oy);
                         else
                         v[0] = 0.0;
                     },
                     [&](reefamr_patch *c, int id, const reefamr_fill &f, const double *w) { Zf(id)(f.di,f.dj) = w[0]; });
        };

        fill_l();

        const int npass = (l<(int)nlevg.size() && nlevg[l]>1) ? RSPASS : 1;
        for(int pass=0; pass<npass; ++pass)
        {
            if(pass>0)
            fill_l();

            for(auto c : rsg)
            {
                if(c->lev!=l)
                continue;

                reefamr_patch *P0 = P[c->g];
                lexer *pp = c->q;
                const matrix2D &M = c->R.e->mat();
                slice &Z = *c->v[z];
                slice &R = *c->v[r];
                sc_level &L = c->pmg->fine();

                for(int ii=EXT; ii<EXT+P0->nx; ++ii)
                for(int jj=EY; jj<EY+P0->ny; ++jj)
                {
                    const long qd = L.idx(ii-EXT,jj-EY,0);
                    L.u[qd] = 0.0;
                    L.f[qd] = 0.0;
                    if(L.act[qd]==0)
                    continue;
                    const int n = c->rowof[lij(pp,ii,jj)];
                    if(n<0)
                    continue;
                    const double az = M.p[n]*Z(ii,jj) + M.n[n]*Z(ii+1,jj) + M.s[n]*Z(ii-1,jj) + M.w[n]*Z(ii,jj+1) + M.e[n]*Z(ii,jj-1);
                    L.f[qd] = R(ii,jj) - az;
                }

                if(sy==0)
                c->xz.solve(L,L.f,L.u);
                else
                c->pmg->vcycle(0,1,1);

                for(int ii=EXT; ii<EXT+P0->nx; ++ii)
                for(int jj=EY; jj<EY+P0->ny; ++jj)
                {
                    const long qd = L.idx(ii-EXT,jj-EY,0);
                    if(L.act[qd])
                    Z(ii,jj) += L.u[qd];
                }
            }
        }
    }
}

void reefamr::rsi_ras(int r, int z)
{
    // one V-cycle per level-0 grid on the residual of its own rows (zero on the cells around it), the result on its
    // own rows: restricted additive Schwarz over the grids with a multigrid subdomain solve
    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;

        const int g = c->g;
        lexer *q = c->q;
        matrix2D &M = const_cast<matrix2D&>(c->R.e->mat());
        slice &R = *c->v[r];
        slice &Z = *c->v[z];
        vec2D &rv = *c->rv;
        const vector<int> &w = c->own;

        std::fill(rv.V.begin(),rv.V.end(),0.0);
        for(size_t k=0; k<w.size(); k+=3)
        rv.V[w[k+2]] = R(w[k],w[k+1]);

        if(g<0)
        c->mg->vcycle(q,rsipgc,*c->z,M,rv,c->rebuild);

        else
        rsi_run(g,[&](){ c->mg->vcycle(q,rsipgc,*c->z,M,rv,c->rebuild); });

        c->rebuild = 0;

        for(size_t k=0; k<w.size(); k+=3)
        Z(w[k],w[k+1]) = (*c->z)(w[k],w[k+1]);
    }
}

// the coarse correction of vector v into the interior of the level-l patches (from the parent of every 2x2 block, on
// its rank, linear interpolation)
void reefamr::rsi_prolong(int l, int v)
{
    auto Vf = [&](int g) -> slice& { for(auto c : rsg) if(c->g==g) return *c->v[v]; return *rsg[0]->v[v]; };

    for(int id : lev[l])
    {
        reefamr_patch *c = P[id];
        slice &x = Vf(id);
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EY; jj<EY+c->ny; ++jj)
        x(ii,jj) = 0.0;
    }

    block_down(l,4,8660+l,
               [&](const reefamr_block &B, int key, double *w)
               {
                   lexer *q = glex(B.g);
                   for(int a=0;a<2;++a)
                   for(int d=0;d<=sy;++d)
                   w[2*a+d] = rsi_interp(q,Vf(B.g),B.ic,B.jc,a==0 ? -1 : 1,(sy==0) ? 0 : (d==0 ? -1 : 1));
               },
               [&](reefamr_patch *q, int n, int kb, const double *w)
               {
                   slice &x = Vf(n);
                   const int nby = q->ny>>sy;
                   const int bi = kb/nby, bj = kb%nby;
                   for(int a=0;a<2;++a)
                   for(int d=0;d<=sy;++d)
                   x(EXT+2*bi+a,EY+ychild(bj)+d) = w[2*a+d];
               });
}

// the patch multigrids of this pass: the rows of every patch (own cells, covered cells included), no coupling across
// the box (the correction is 0 there)
void reefamr::rsi_pmg()
{
    for(auto c : rsg)
    {
        if(c->lev==0)
        continue;

        reefamr_patch *P0 = P[c->g];
        lexer *pp = c->q;
        const matrix2D &M = c->R.e->mat();

        if(c->pmg==nullptr)
        {
            c->pmg = new reefmg_core;
            c->pmg->set_precision(64);
            if(!c->pmg->setup(MPI_COMM_SELF,P0->nx,P0->ny,1,P0->nx,P0->ny,0,pp->DXN+marge+EXT-1,pp->DYN+marge+EY-1))
            cout<<"REEFAMR semi-implicit: patch multigrid "<<c->pmg->err()<<endl;
        }

        sc_level &L = c->pmg->fine();
        std::fill(L.p.begin(),L.p.end(),0.0);
        std::fill(L.n.begin(),L.n.end(),0.0); std::fill(L.s.begin(),L.s.end(),0.0);
        std::fill(L.w.begin(),L.w.end(),0.0); std::fill(L.e.begin(),L.e.end(),0.0);
        std::fill(L.t.begin(),L.t.end(),0.0); std::fill(L.b.begin(),L.b.end(),0.0);
        std::fill(L.u.begin(),L.u.end(),0.0); std::fill(L.f.begin(),L.f.end(),0.0);
        std::fill(L.act.begin(),L.act.end(),0);

        // fill the inactive cells (solids) with identity rows
        for(int ii=EXT; ii<EXT+P0->nx; ++ii)
        for(int jj=EY; jj<EY+P0->ny; ++jj)
        L.p[L.idx(ii-EXT,jj-EY,0)] = 1.0;

        double ldx = 0.0;
        auto inner = [&](int ii, int jj) { return ii>=EXT && ii<EXT+P0->nx && jj>=EY && jj<EY+P0->ny; };

        for(size_t q=0; q<c->own.size(); q+=3)
        {
            const int ii=c->own[q], jj=c->own[q+1], n=c->own[q+2];
            const long qd = L.idx(ii-EXT,jj-EY,0);

            double cf[4] = {M.n[n],M.s[n],M.w[n],M.e[n]};
            const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
            for(int d=0; d<4; ++d)
            {
                const int a = ii+di[d], b = jj+dj[d];
                if(!inner(a,b) || pp->flagslice4[lij(pp,a,b)]<0)
                cf[d] = 0.0;
            }

            const double mp = M.p[n];
            const double react = mp+M.n[n]+M.s[n]+M.w[n]+M.e[n];
            const double diff = fabs(cf[0])+fabs(cf[1])+fabs(cf[2])+fabs(cf[3]);
            if(react>1.0e-12*fabs(mp))
            ldx = MAX(ldx,sqrt(diff/(4.0*react)));

            L.p[qd]=mp; L.n[qd]=cf[0]; L.s[qd]=cf[1]; L.w[qd]=cf[2]; L.e[qd]=cf[3];
            L.act[qd]=1;
        }

        const int need = 2 + (int)ceil(log2(ldx>1.0 ? ldx : 1.0));
        c->pmg->set_active_levels(need);

        if(sy==0)
        c->xz.factor(L,64);
        else
        c->pmg->coarsen();
    }
}

// ---- G 57: the coarse space of the eta solve across the level-0 grids (native grid and blocks): aggregates of
// G 57 x G 57 own columns of every grid, the Galerkin operator of the piecewise constant prolongation (the rows
// summed over the aggregates), factorised on every rank
void reefamr::rsi_cs_build()
{
    rsics_n = 0;
    const int csw = p0->G57;
    if(csw<=0 || nblk==0)
    return;

    if(rsics==nullptr)
    rsics = new rsi_csi;

    // the aggregates of the own rows, numbered over the ranks
    int nloc = 0;
    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;
        lexer *q = c->q;
        slice &AG = *c->v[VAG];
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        AG(ii,jj) = -1.0;

        std::map<long,int> ids;
        const vector<int> &w = c->own;
        for(size_t k=0; k<w.size(); k+=3)
        {
            const long key = long(w[k]/csw)*100000L + long(w[k+1]/csw);
            auto it = ids.find(key);
            if(it==ids.end())
            it = ids.emplace(key,nloc++).first;
            AG(w[k],w[k+1]) = it->second;
        }
    }

    int off = 0;
    MPI_Exscan(&nloc,&off,1,MPI_INT,MPI_SUM,reef3d_comm());
    if(p0->mpirank==0)
    off = 0;
    int ntot = 0;
    MPI_Allreduce(&nloc,&ntot,1,MPI_INT,MPI_SUM,reef3d_comm());

    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;
        slice &AG = *c->v[VAG];
        const vector<int> &w = c->own;
        for(size_t k=0; k<w.size(); k+=3)
        AG(w[k],w[k+1]) += off;
    }

    // the aggregates of the cells around the blocks, of the ring and of the halo
    rsi_fillv(VAG,false);

    // the coarse operator
    std::map<std::pair<int,int>,double> T;
    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;
        const matrix2D &M = c->R.e->mat();
        slice &AG = *c->v[VAG];
        const vector<int> &w = c->own;
        for(size_t k=0; k<w.size(); k+=3)
        {
            const int ii=w[k], jj=w[k+1], n=w[k+2];
            const int I = int(AG(ii,jj));
            T[{I,I}] += M.p[n];
            const double cf[4] = {M.n[n],M.s[n],M.w[n],M.e[n]};
            const double an[4] = {AG(ii+1,jj),AG(ii-1,jj),AG(ii,jj+1),AG(ii,jj-1)};
            for(int m=0; m<4; ++m)
            if(cf[m]!=0.0 && an[m]>=0.0 && an[m]==an[m])
            T[{I,int(an[m])}] += cf[m];
        }
    }

    vector<int> ti, tj;
    vector<double> tv;
    for(auto &e : T)
    {
        ti.push_back(e.first.first);
        tj.push_back(e.first.second);
        tv.push_back(e.second);
    }
    const int np = p0->mpi_size;
    int nt = int(tv.size());
    vector<int> cnt(np), dsp(np,0);
    MPI_Allgather(&nt,1,MPI_INT,&cnt[0],1,MPI_INT,reef3d_comm());
    for(int r=1; r<np; ++r)
    dsp[r] = dsp[r-1] + cnt[r-1];
    const int nall = dsp[np-1] + cnt[np-1];
    vector<int> ai(MAX(nall,1)), aj(MAX(nall,1));
    vector<double> av(MAX(nall,1));
    MPI_Allgatherv(ti.data(),nt,MPI_INT,&ai[0],&cnt[0],&dsp[0],MPI_INT,reef3d_comm());
    MPI_Allgatherv(tj.data(),nt,MPI_INT,&aj[0],&cnt[0],&dsp[0],MPI_INT,reef3d_comm());
    MPI_Allgatherv(tv.data(),nt,MPI_DOUBLE,&av[0],&cnt[0],&dsp[0],MPI_DOUBLE,reef3d_comm());

    if(ntot==0)
    return;

    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(nall);
    for(int k=0; k<nall; ++k)
    trip.emplace_back(ai[k],aj[k],av[k]);
    Eigen::SparseMatrix<double> Ac(ntot,ntot);
    Ac.setFromTriplets(trip.begin(),trip.end());
    Ac.makeCompressed();
    rsics->lu.compute(Ac);
    if(rsics->lu.info()!=Eigen::Success)
    {
        if(p0->mpirank==0)
        cout<<"semi-implicit multi-block: the coarse space of the eta solve (G 57) is singular, off for this solve"<<endl;
        return;
    }
    rsics->rl.assign(ntot,0.0);
    rsics->rc.assign(ntot,0.0);
    rsics_n = ntot;
}

// z = P Ac^-1 R r on the own rows, 0 on all other cells of the level-0 grids
void reefamr::rsi_cs_correct(int r, int z)
{
    std::fill(rsics->rl.begin(),rsics->rl.end(),0.0);

    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;
        slice &R = *c->v[r];
        slice &AG = *c->v[VAG];
        const vector<int> &w = c->own;
        for(size_t k=0; k<w.size(); k+=3)
        rsics->rl[int(AG(w[k],w[k+1]))] += R(w[k],w[k+1]);
    }

    MPI_Allreduce(rsics->rl.data(),rsics->rc.data(),rsics_n,MPI_DOUBLE,MPI_SUM,reef3d_comm());

    Eigen::Map<Eigen::VectorXd> rc(rsics->rc.data(),rsics_n);
    Eigen::VectorXd zc = rsics->lu.solve(rc);

    for(auto c : rsg)
    {
        if(c->lev>0)
        continue;
        lexer *q = c->q;
        slice &Z = *c->v[z];
        slice &AG = *c->v[VAG];
        for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
        for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
        Z(ii,jj) = 0.0;
        const vector<int> &w = c->own;
        for(size_t k=0; k<w.size(); k+=3)
        Z(w[k],w[k+1]) = zc[int(AG(w[k],w[k+1]))];
    }
}

// ------------------------------------------------------------------------------------------- the solve
// the coefficients a_f of the coarse faces next to the patches: the mean of the two fine ones (finest level first, so
// that a patch with children sends the values it took)
void reefamr::rsi_couple()
{
    rsm_begin(1);

    for(int l=maxlev; l>=0; --l)
    {
        // the coarse faces of the grids of level l take the values of their children
        for(auto c : rsg)
        {
            if(c->lev!=l)
            continue;

            slice &ax = *c->R.ax;
            slice &ay = *c->R.ay;

            rsm_take(c->g,[&](int dir, int fi, int fj, const double *v, int sgn)
            {
                if(dir==1)
                ax(fi,fj) = v[0];
                else
                ay(fi,fj) = v[0];
            });
        }

        if(l==0)
        break;

        // the box faces of the patches of level l
        for(auto c : rsg)
        {
            if(c->lev!=l)
            continue;

            reefamr_patch *P0 = P[c->g];
            slice &ax = *c->R.ax;
            slice &ay = *c->R.ay;
            const int il = EXT-1, ih = EXT+P0->nx-1, jl = EY-1, jh = EY+P0->ny-1;

            for(int r=0; r<P0->ny; ++r)
            {
                *rsm_slot(c->g,0,r) = ax(il,EY+r);
                *rsm_slot(c->g,1,r) = ax(ih,EY+r);
            }

            if(sy==1)
            for(int r=0; r<P0->nx; ++r)
            {
                *rsm_slot(c->g,2,r) = ay(EXT+r,jl);
                *rsm_slot(c->g,3,r) = ay(EXT+r,jh);
            }
        }

        rsm_send(l,8580);
    }
}

void reefamr::rsi_solve(lexer *p, ghostcell *pgc, int itmax, rsi_result &res)
{
    rsipgc = pgc;
    rsi_grids();

    const double t0 = pgc->timer();
    int it_done=0, iter=0;

    for(auto c : rsg)
    c->R.e->wd_start(c->q);

    // a coarse face next to a finer patch takes the mean coefficient a_f of its fine faces: the operator replaces its
    // flux by the fine one anyway (rsi_apply), but the coefficient decides whether the coarse cell is coupled (the
    // wetting and drying pattern of si_eta2D::wd_update, the isolated rows of si_eta2D::assemble). With the coarse
    // coefficient a coarse cell whose own faces are closed would keep its row without the fine flux, while the face
    // update gives it the fine discharge: its water column would go below zero and be clipped (volume gained)
    if(rsi_lev)
    rsi_couple();

    rsi_space sp{this};

    auto etaF = [&](int g) -> slice& { for(auto c : rsg) if(c->g==g) return *c->R.etaI; return *rsg[0]->R.etaI; };

    for(int it=0; it<itmax; ++it)
    {
        // the rows of every grid, the right-hand side and the initial guess
        for(auto c : rsg)
        {
            lexer *q = c->q;
            si_eta2D *e = c->R.e;

            e->wd_rows(q,*c->R.rhsW,*c->R.depth,*c->R.bdiag);
            e->build(q,*c->R.ax,*c->R.ay,e->wd_rhs(),*c->R.etaI,e->wd_react());

            const vec2D &b = e->rhsv();
            const vector<int> &w = c->own;
            for(size_t k=0; k<w.size(); k+=3)
            {
            (*c->v[VB])(w[k],w[k+1]) = b.V[w[k+2]];
            (*c->v[VX])(w[k],w[k+1]) = (*c->R.etaI)(w[k],w[k+1]);
            }
        }

        // the matrices changed: the preconditioner coarsens again at its first use (keeping the
        // hierarchy over the passes of the wetting and drying iteration doubles the iterations)
        for(auto c : rsg)
        c->rebuild = 1;

        // G 57: the coarse operator of this pass; the patch multigrids
        rsi_cs_build();
        if(rsi_lev)
        rsi_pmg();

        rsi_apply(VX,VVV);

        double bn,rn;
        int restarts=0;
        iter = reefamr_bicgstab(sp,p->N44,p->N46,bn,rn,nullptr,nullptr,&restarts);

        if(iter>=p->N46 && p->mpirank==0)
        cout<<"semi-implicit, all grids: eta solve at the iteration limit N 46 = "<<p->N46<<", residual "<<(bn>0.0?rn/bn:rn)<<endl;

        // eta^{n+1} of the leaf cells, then of the covered cells and around the grids
        for(auto c : rsg)
        {
            const vector<int> &w = c->leaf;
            for(size_t k=0; k<w.size(); k+=3)
            (*c->R.etaI)(w[k],w[k+1]) = (*c->v[VX])(w[k],w[k+1]);
        }

        rsi_fill(7630,etaF);
        rsi_bound();

        ++it_done;

        if(itmax==1)
        break;

        // new wet pattern
        int change=0;
        for(auto c : rsg)
        change += c->R.e->wd_update(c->q,*c->R.ax,*c->R.ay,*c->R.bdiag,*c->R.etaI,*c->R.depth);

        int all=0;
        MPI_Allreduce(&change,&all,1,MPI_INT,MPI_SUM,pgc->mpi_comm);

        if(all==0)
        break;
    }

    // dry cells
    int dry=0, alldry=0;
    for(auto c : rsg)
    {
        si_eta2D *e = c->R.e;
        const vector<int> &w = c->leaf;
        for(size_t k=0; k<w.size(); k+=3)
        if(e->wd_pattern(w[k],w[k+1])<0.5)
        ++dry;
    }
    MPI_Allreduce(&dry,&alldry,1,MPI_INT,MPI_SUM,pgc->mpi_comm);

    res.iter = iter;
    res.wd_iter = it_done;
    res.dry = alldry;
    res.time = pgc->timer()-t0;
}

// ------------------------------------------------------------------------------------------- face matching
int reefamr::rsm_count(int id, int side)
{
    reefamr_patch *c = P[id];
    if(side<2)
    return c->ny;
    return (sy==1) ? c->nx : 0;
}

void reefamr::rsm_begin(int nv)
{
    rsm_nv = nv;
    rsm_rec.resize(P.size());

    for(size_t id=0; id<P.size(); ++id)
    {
        if(P[id]->lev<1)
        continue;
        const size_t n = (size_t)(2*rsm_count(id,0) + 2*rsm_count(id,2))*nv;
        rsm_rec[id].assign(n,0.0);
    }

    rsm_val.resize(rmatch.size());
    for(size_t t=0; t<rmatch.size(); ++t)
    rsm_val[t].assign(rmatch[t].size()*nv,0.0);
}

double* reefamr::rsm_slot(int id, int side, int r)
{
    const int ny = rsm_count(id,0), nx = rsm_count(id,2);
    const int o[4] = {0, ny, 2*ny, 2*ny+nx};
    return &rsm_rec[id][(size_t)(o[side]+r)*rsm_nv];
}

void reefamr::rsm_face(int id, int side, int r, int &dir, int &fi, int &fj)
{
    reefamr_patch *c = P[id];
    if(side==0) { dir=1; fi=EXT-1;       fj=EY+r; }
    if(side==1) { dir=1; fi=EXT+c->nx-1; fj=EY+r; }
    if(side==2) { dir=2; fi=EXT+r;       fj=EY-1; }
    if(side==3) { dir=2; fi=EXT+r;       fj=EY+c->ny-1; }
}

void reefamr::rsm_send(int l, int tag)
{
    const int nv = rsm_nv;

    face_run_at(l,nv,tag+l,
                [&](reefamr_patch *q, int side, int r, double *sb)
                {
                    int id = -1;
                    for(int n : lev[l])
                    if(P[n]==q)
                    {
                        id = n;
                        break;
                    }
                    const double *a = rsm_slot(id,side,r);
                    const double *b = rsm_slot(id,side,r+fr1(side));
                    for(int v=0; v<nv; ++v)
                    sb[v] = 0.5*(a[v]+b[v]);
                },
                [&](int tgt, int idx, const double *w)
                {
                    for(int v=0; v<nv; ++v)
                    rsm_val[tgt][(size_t)idx*nv+v] = w[v];
                });
}

void reefamr::rsm_take(int g, const std::function<void(int, int, int, const double*, int)> &f)
{
    const int t = g+1;
    if(t>=(int)match.size())
    return;

    const int nv = rsm_nv;
    vector<double> m(nv);

    for(auto &e : match[t])
    {
        const double *a = rsm_slot(e.child,e.side,e.r);
        const double *b = rsm_slot(e.child,e.side,e.r+fr1(e.side));
        for(int v=0; v<nv; ++v)
        m[v] = 0.5*(a[v]+b[v]);
        f(e.dir+1,e.fi,e.fj,m.data(),e.sgn);
    }

    for(size_t k=0; k<rmatch[t].size(); ++k)
    {
        reefamr_match &e = rmatch[t][k];
        f(e.dir+1,e.fi,e.fj,&rsm_val[t][k*nv],e.sgn);
    }
}
