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

#include"fnpf_amr.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice4.h"
#include"fnpf_fsfbc.h"
#include"fnpf_sigma.h"
#include"fnpf_fsf_update.h"
#include"fnpf_bed_update.h"
#include"fnpf_laplace_cds2.h"
#include"reefmg_core.h"
#include"fnpf_amr_fill.h"
#include"fnpf_body.h"
#include<cmath>
#include<mpi.h>
#include<iomanip>
#include<cstdio>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
inline int lij(const lexer *q, int ii, int jj) { return (ii-q->imin)*q->jmax + (jj-q->jmin); }

typedef reefamr_comms_off comms_off;

// the Laplace solver seen by fnpf_RK3 with mesh refinement: level 0 is assembled by its own
// solver (assemble_only), then fnpf_amr solves all grids together
class fnpf_laplace_amr : public fnpf_laplace
{
public:
    fnpf_laplace_amr(fnpf_amr *a, fnpf_laplace *outer, fnpf_laplace_cds2 *plap0) : a(a), outer(outer), plap0(plap0) {}

    void start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, fnpf_fsf *pf, double *f, slice &Fifsf) override
    {
        if(!a->active())
        {
            outer->start(p,c,pgc,psolv,pf,f,Fifsf);
            return;
        }

        a->lap_solve(p,c,pgc,psolv,pf,f,Fifsf);
    }

private:
    fnpf_amr *a;
    fnpf_laplace *outer;
    fnpf_laplace_cds2 *plap0;
};
}

fnpf_amr::fnpf_amr(lexer *p, fdm_fnpf *c, ghostcell *pgc) : reefamr(p,pgc)
{
    c0 = c;

    reefamr_param q;
    q.name = "FNPF AMR";
    q.maxlev = p->A270;
    q.regrid = 0;
    q.nbuf = MAX(p->A272,0);
    q.tile = MAX(p->A275,4);
    q.tile += q.tile%2;
    q.keep = 0;

    // cells computed beyond the patch box: the free-surface derivatives (WENO5, CDS4) reach
    // three cells; the patch arrays reach EXT + margin fine cells beyond the patch, which the
    // coarser level has to cover
    q.ext = 3;
    q.nest = MAX(3,(q.ext+p->margin+1)/2);

    // vertical refinement: every refined level doubles the sigma layers (A 281 1)
    vref = (p->A281==1) ? 2 : 1;
    q.vref.assign(q.maxlev+1,vref);

    for(int k=0; k<p->A276; ++k)
    {
        q.rbox.push_back(p->A276_xs[k]); q.rbox.push_back(p->A276_xe[k]);
        q.rbox.push_back(p->A276_ys[k]); q.rbox.push_back(p->A276_ye[k]);
    }
    for(int k=0; k<p->A277; ++k)
    {
        q.fbox.push_back(p->A277_xs[k]); q.fbox.push_back(p->A277_xe[k]);
        q.fbox.push_back(p->A277_ys[k]); q.fbox.push_back(p->A277_ye[k]);
    }

    // no refinement in the relaxation zones of the wave generation (B 98 2) and the
    // numerical beach (B 99 1, 2), measured from the ends of the domain in x
    const double big = 1.0e20;
    if(p->B98==2 && p->B96_1>0.0)
    {
        q.fbox.push_back(-big); q.fbox.push_back(p->global_xmin+p->B96_1);
        q.fbox.push_back(-big); q.fbox.push_back(big);
    }
    if((p->B99==1 || p->B99==2) && p->B96_2>0.0)
    {
        q.fbox.push_back(p->global_xmax-p->B96_2); q.fbox.push_back(big);
        q.fbox.push_back(-big); q.fbox.push_back(big);
    }
    q.ioband = 4;

    // refinement around the resolved body (X 10 1): margin A 278 around the wetted hull
    q.zones = (p->X10==1 && p->A278>0);
    q.zr = p->A278_r;
    q.zL = p->A279_L;
    q.za = p->A279_a;

    configure(q);

    gcval_eta = 55;
    gcval_fifsf = 60;

    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;
    stg = 2;

    printtime_amr = 0.0;
    printcount_amr = 0;
    lap_it_total = lap_solves = 0;
    lap_it_last = 0;
    lap_res_last = 0.0;
    for(int k=0; k<6; ++k)
    tm[k] = 0.0;
}

fnpf_amr::~fnpf_amr()
{
    free_patches();

    for(auto v : kv0)
    delete [] v;
    delete mg0;
}

fnpf_laplace* fnpf_amr::laplace(fnpf_laplace *outer, fnpf_laplace_cds2 *plap0)
{
    if(maxlev<1)
    return outer;

    fnpf_amr::plap0 = plap0;

    return new fnpf_laplace_amr(this,outer,plap0);
}

void fnpf_amr::attach_body(fnpf_body *b)
{
    body = (b!=nullptr && b->present()) ? b : nullptr;
}

void fnpf_amr::zone_bodies(vector<sixdof_obj*> &obj)
{
    if(body!=nullptr)
    body->amr_bodies(obj);
}

fnpf_fsf* fnpf_amr::patch_fsf(int n)
{
    return FP(n)->pf;
}

slice& fnpf_amr::patch_tendency(int n, int m)
{
    return (m==0) ? static_cast<slice&>(*FP(n)->ek) : static_cast<slice&>(*FP(n)->fk);
}

// the finest local grid whose interior holds (x,y)
int fnpf_amr::finest_at(double x, double y)
{
    int g=-1, l=0;
    for(int n=0; n<(int)P.size(); ++n)
    {
        reefamr_patch *c = P[n];
        lexer *pp = c->pp;
        if(c->lev<=l)
        continue;
        if(x>=pp->XN[EXT+marge] && x<pp->XN[EXT+c->nx+marge] && y>=pp->YN[EXT+marge] && y<pp->YN[EXT+c->ny+marge])
        {
            g=n;
            l=c->lev;
        }
    }
    return g;
}

// --------------------------------------------------------------------- grids
int fnpf_amr::fidx(lexer *q, int ii, int jj, int kk) const
{
    return (ii-q->imin)*q->jmax*q->kmaxF + (jj-q->jmin)*q->kmaxF + kk - q->kmin;
}

int fnpf_amr::gknoz(int g)
{
    return glex(g)->knoz;
}

slice& fnpf_amr::sval(int g, int s, int m)
{
    if(g<0)
    return (m==0) ? *Se0 : *Sf0;

    fnpf_amr_patch *c = FP(g);
    if(s==0) return (m==0) ? *c->erk1 : *c->frk1;
    if(s==1) return (m==0) ? *c->erk2 : *c->frk2;
    return (m==0) ? c->c->eta : c->c->Fifsf;
}

// --------------------------------------------------------------------- setup
void fnpf_amr::ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    int ok=1;
    if(p->A310!=3) ok=0;
    if(p->A343!=0) ok=0;
    if(p->A350!=0) ok=0;
    if(p->X10!=0 && p->X10!=1) ok=0;
    if(p->A324!=0) ok=0;
    if(p->A328!=0) ok=0;
    if(p->A380!=0) ok=0;
    if(p->A311==3) ok=0;
    if(p->j_dir!=1) ok=0;

    if(ok==0)
    {
        if(p->mpirank==0)
        cout<<"FNPF AMR (A 270): only for A 310 3, A 311 other than 3, A 343 0, A 350 0, X 10 0/1, A 324 0, A 328 0, A 380 0 and 3D grids -- refinement switched off"<<endl;
        maxlev=0;
        return;
    }

    setup(p,pgc);

    mkdir("./REEF3D_FNPF_AMR",0777);

    if(p->mpirank==0)
    {
        logout.open("./REEF3D_FNPF_AMR/REEF3D_FNPF_AMR_log.dat");
        logout<<"# count \t simtime \t dt \t patches \t cells \t Laplace iterations \t residual"<<endl;
    }

    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;
    stg = 2;

    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);

    if(body!=nullptr)
    body->amr_grids(p,pgc);

    if(p->mpirank==0)
    {
        cout<<"FNPF AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined columns";
        if(vref==2)
        cout<<", sigma layers doubled on every level (A 281)";
        cout<<endl;
    }
}

// --------------------------------------------------------------------- patches
reefamr_patch* fnpf_amr::patch_new()
{
    return new fnpf_amr_patch;
}

// sigma grid and 3D flags of the patch, as driver::makegrid_sigma: a fluid column is fluid
// from the bed node k=0 to the free-surface node k=knoz
void fnpf_amr::build_lexer3D(fnpf_amr_patch &c)
{
    lexer *pp = c.pp;
    const int n3 = pp->imax*pp->jmax*pp->kmax;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);

    pp->Iarray(pp->flag4,n3);
    pp->Iarray(pp->flag7,n7);
    for(int q=0; q<n7; ++q)
    pp->flag7[q] = -10;

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        const bool fl = pp->flagslice4[lij(pp,ii,jj)]>0;

        for(int kk=pp->kmin; kk<pp->kmin+pp->kmax; ++kk)
        pp->flag4[(ii-pp->imin)*pp->jmax*pp->kmax + (jj-pp->jmin)*pp->kmax + kk-pp->kmin] = (fl && kk>=0 && kk<pp->knoz) ? 1 : -10;

        if(fl)
        for(int kk=0; kk<=pp->knoz; ++kk)
        pp->flag7[fidx(pp,ii,jj,kk)] = 1;
    }

    pp->Darray(pp->ZSN,pp->imax*pp->jmax*(pp->kmax+1));
    pp->Darray(pp->ZSP,n3);
    pp->Darray(pp->bed,pp->imax*pp->jmax);

    // the free surface spans the whole patch: no inflow clamp of velcalc_sig, no domain ends
    pp->origin_i = pp->origin_j = 100000000;
}

void fnpf_amr::free_lexer3D(fnpf_amr_patch &c)
{
    lexer *pp = c.pp;
    const int n3 = pp->imax*pp->jmax*pp->kmax;
    const int n7 = pp->imax*pp->jmax*(pp->kmax+2);

    pp->del_Iarray(pp->flag4,n3);
    pp->del_Iarray(pp->flag7,n7);
    pp->del_Darray(pp->ZSN,pp->imax*pp->jmax*(pp->kmax+1));
    pp->del_Darray(pp->ZSP,n3);
    pp->del_Darray(pp->bed,pp->imax*pp->jmax);
    pp->del_Darray(pp->sig,n7);
    pp->del_Darray(pp->sigx,n7);
    pp->del_Darray(pp->sigy,n7);
    pp->del_Darray(pp->sigz,pp->imax*pp->jmax);
    pp->del_Darray(pp->sigxx,n7);
}

// columns of the patch with a solid neighbour in x or y (wall ghost nodes of Fi, as gc_fivec)
void fnpf_amr::build_walls(fnpf_amr_patch &c)
{
    lexer *pp = c.pp;
    c.wall.clear();

    for(int ii=0; ii<pp->knox; ++ii)
    for(int jj=0; jj<pp->knoy; ++jj)
    {
        if(pp->flagslice4[lij(pp,ii,jj)]<0)
        continue;

        if(pp->flagslice4[lij(pp,ii-1,jj)]<0 || pp->flagslice4[lij(pp,ii+1,jj)]<0
        || pp->flagslice4[lij(pp,ii,jj-1)]<0 || pp->flagslice4[lij(pp,ii,jj+1)]<0)
        {
            c.wall.push_back(ii);
            c.wall.push_back(jj);
        }
    }
}

// wall ghost nodes of a Fi-layout array: mirrored, as ghostcell::fivec (no in/outflow on a patch)
void fnpf_amr::walls_fi(fnpf_amr_patch &c, double *f)
{
    lexer *pp = c.pp;
    const int sI = pp->jmax*pp->kmaxF;
    const int sJ = pp->kmaxF;
    const int *fl = pp->flag7;

    for(size_t n=0; n<c.wall.size(); n+=2)
    for(int kk=0; kk<=pp->knoz; ++kk)
    {
        const int q = fidx(pp,c.wall[n],c.wall[n+1],kk);

        if(fl[q-sI]<0) { f[q-sI] = f[q]; f[q-2*sI] = f[q]; f[q-3*sI] = f[q]; }
        if(fl[q+sI]<0) { f[q+sI] = f[q]; f[q+2*sI] = f[q]; f[q+3*sI] = f[q]; }
        if(fl[q-sJ]<0) { f[q-sJ] = f[q]; f[q-2*sJ] = f[q]; f[q-3*sJ] = f[q]; }
        if(fl[q+sJ]<0) { f[q+sJ] = f[q]; f[q+2*sJ] = f[q]; f[q+3*sJ] = f[q]; }
    }
}

void fnpf_amr::walls_sl(fnpf_amr_patch &c, slice &f, int gcv)
{
    comms_off guard(pgc0);
    pgc0->gcsl_start4(c.pp,f,gcv);
}

// FNPF objects of a new patch (comms off); bed and state are set by regrid_state
void fnpf_amr::patch_objects(reefamr_patch *q, ghostcell *pgc)
{
    fnpf_amr_patch *c = FP(q);

    build_bc2D(pgc,*c);
    build_lexer3D(*c);

    lexer *pp = c->pp;

    c->c = new fdm_fnpf(pp);
    c->pf = new fnpf_fsfbc(pp,c->c,pgc);
    c->psig = new fnpf_sigma(pp,c->c,pgc);
    c->pfu = new fnpf_fsf_update(pp,c->c,pgc);
    c->pbu = new fnpf_bed_update(pp);
    c->plap = new fnpf_laplace_cds2(pp);
    c->plap->assemble_only = true;

    c->erk1 = new slice4(pp); c->erk2 = new slice4(pp);
    c->frk1 = new slice4(pp); c->frk2 = new slice4(pp);
    c->ek = new slice4(pp); c->fk = new slice4(pp);

    fdm_fnpf *cc = c->c;
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
        cc->bc(ii,jj) = 0;
        cc->etaloc(ii,jj) = pp->knoz;
        pp->wet[lij(pp,ii,jj)] = 1;
    }

    build_walls(*c);
}

void fnpf_amr::patch_delete(reefamr_patch *q)
{
    fnpf_amr_patch *c = FP(q);

    for(auto v : c->kv)
    delete [] v;
    delete c->mg;

    delete c->erk1; delete c->erk2; delete c->frk1; delete c->frk2; delete c->ek; delete c->fk;
    delete c->plap;
    delete c->pbu;
    delete c->pfu;
    delete c->psig;
    delete c->pf;
    delete c->c;

    free_lexer3D(*c);
}

void fnpf_amr::tag(int l, vector<unsigned char> &M)
{
    // static refinement: the boxes (A 276) are marked by the core
}

// --------------------------------------------------------------------- interpolation
// interpolation from coarse cell (ic,jc) to the centre of its child (ox,oy), a quarter coarse
// cell away in x and y: weights w[(di+2)*5+(dj+2)] of the coarse cells (ic+di,jc+dj).
// pord 4: bicubic (tensor product, the 4 cells around the child centre in each direction),
// where all 16 cells are fluid; otherwise biquadratic on the 3x3 cells, slopes and
// curvatures switched off next to solids
void fnpf_amr::pweights(lexer *q, int ic, int jc, int ox, int oy, double *w)
{
    auto fl = [&](int a, int b) { return q->flagslice4[lij(q,a,b)]>0; };

    for(int m=0; m<25; ++m)
    w[m] = 0.0;

    if(pord>=4)
    {
        bool all=true;
        const int a0 = (ox>0) ? -1 : -2, b0 = (oy>0) ? -1 : -2;
        for(int a=a0; a<a0+4 && all; ++a)
        for(int b=b0; b<b0+4 && all; ++b)
        if(!fl(ic+a,jc+b))
        all=false;

        if(all)
        {
            // Lagrange weights at +1/4 on the nodes -1,0,1,2 (mirrored for -1/4)
            const double cp[4] = {-0.0546875, 0.8203125, 0.2734375, -0.0390625};
            double wx[4], wy[4];
            for(int m=0; m<4; ++m)
            {
                wx[m] = (ox>0) ? cp[m] : cp[3-m];
                wy[m] = (oy>0) ? cp[m] : cp[3-m];
            }
            for(int a=0; a<4; ++a)
            for(int b=0; b<4; ++b)
            w[(a0+a+2)*5+(b0+b+2)] = wx[a]*wy[b];
            return;
        }
    }

    const bool ex = fl(ic+1,jc) && fl(ic-1,jc);
    const bool ey = fl(ic,jc+1) && fl(ic,jc-1);
    const bool exy = ex && ey && fl(ic+1,jc+1) && fl(ic-1,jc+1) && fl(ic+1,jc-1) && fl(ic-1,jc-1);

    const double x = 0.25*ox, y = 0.25*oy;
    auto W = [&](int di, int dj) -> double& { return w[(di+2)*5+(dj+2)]; };

    W(0,0) = 1.0;

    if(ex)
    {
        W(1,0) += 0.5*x + 0.5*x*x;
        W(-1,0) += -0.5*x + 0.5*x*x;
        W(0,0) -= x*x;
    }
    if(ey)
    {
        W(0,1) += 0.5*y + 0.5*y*y;
        W(0,-1) += -0.5*y + 0.5*y*y;
        W(0,0) -= y*y;
    }
    if(exy)
    {
        const double c = 0.25*x*y;
        W(1,1) += c; W(-1,1) -= c; W(1,-1) -= c; W(-1,-1) += c;
    }
}

// slice value at the child centre
double fnpf_amr::pq(slice &f, lexer *q, int ic, int jc, int ox, int oy)
{
    double w[25];
    pweights(q,ic,jc,ox,oy,w);

    double r = 0.0;
    for(int di=-2; di<=2; ++di)
    for(int dj=-2; dj<=2; ++dj)
    {
        const double c = w[(di+2)*5+(dj+2)];
        if(c!=0.0)
        r += c*f(ic+di,jc+dj);
    }
    return r;
}

// column of a child of coarse cell (ic,jc) of grid g: nodes 0..knf of the fine column;
// horizontally biquadratic on the coarse nodes, vertically linear in sigma (the fine nodes
// subdivide the coarse ones evenly)
void fnpf_amr::pcol(int g, int ic, int jc, int ox, int oy, const double *src, int knf, double *v)
{
    lexer *q = glex(g);
    const int knc = q->knoz;
    const int fz = knf/knc;
    const int sI = q->jmax*q->kmaxF;
    const int sJ = q->kmaxF;

    double w[25];
    pweights(q,ic,jc,ox,oy,w);

    const double *s[25];
    int nw=0;
    double ww[25];
    const int n0 = fidx(q,ic,jc,0);
    for(int di=-2; di<=2; ++di)
    for(int dj=-2; dj<=2; ++dj)
    {
        const double c = w[(di+2)*5+(dj+2)];
        if(c==0.0)
        continue;
        s[nw] = src + n0 + di*sI + dj*sJ;
        ww[nw] = c;
        ++nw;
    }

    auto node = [&](int K)
    {
        double r = 0.0;
        for(int m=0; m<nw; ++m)
        r += ww[m]*s[m][K];
        return r;
    };

    if(fz==1)
    {
        for(int kk=0; kk<=knf; ++kk)
        v[kk] = node(kk);
        return;
    }

    double lo = node(0);
    for(int K=0; K<knc; ++K)
    {
        const double hi = node(K+1);
        for(int r=0; r<fz; ++r)
        v[K*fz+r] = lo + (hi-lo)*(double(r)/double(fz));
        lo = hi;
    }
    v[knf] = lo;
}

// --------------------------------------------------------------------- regridding
void fnpf_amr::regrid_prepare(ghostcell *pgc)
{
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;

    auto sel = [&](int g, int m) -> slice& { return sval(g,2,m); };
    for(int l=1; l<=maxlev; ++l)
    fill_sl(l,2,7500+l,sel);
}

void fnpf_amr::regrid_static(ghostcell *pgc)
{
}

// state of the new patches, coarse to fine: bed, depth and its derivatives, eta and Fifsf,
// sigma grid, Fi and Fz from the parent level
void fnpf_amr::regrid_state(ghostcell *pgc, vector<reefamr_patch*> &oldP)
{
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;

    auto sbed = [&](int g, int m) -> slice& { return gfd(g)->bed; };
    auto sfz  = [&](int g, int m) -> slice& { return gfd(g)->Fz; };
    auto sst  = [&](int g, int m) -> slice& { return sval(g,2,m); };
    auto sfi  = [&](int g) -> double* { return gfd(g)->Fi; };

    for(int l=1; l<=maxlev; ++l)
    {
        // bed
        for(int id : lev[l])
        if(P[id]->fresh)
        prolong_interior_sl(*FP(id),1,sbed);
        fill_sl(l,1,7510+l,sbed);

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            lexer *p = c->pp;
            fdm_fnpf *cc = c->c;

            walls_sl(*c,cc->bed,50);

            comms_off guard(pgc);

            SLICEBASELOOP
            {
            p->bed[IJ] = cc->bed(i,j);
            cc->depth(i,j) = p->wd - cc->bed(i,j);
            }
            pgc->gcsl_start4(p,cc->depth,50);

            for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
            for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
            p->bed[lij(p,ii,jj)] = cc->bed(ii,jj);
        }

        // surface
        for(int id : lev[l])
        if(P[id]->fresh)
        {
            prolong_interior_sl(*FP(id),2,sst);
            prolong_interior_sl(*FP(id),1,sfz);
        }
        fill_sl(l,2,7520+l,sst);
        fill_sl(l,1,7530+l,sfz);

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            lexer *p = c->pp;
            fdm_fnpf *cc = c->c;

            walls_sl(*c,cc->eta,gcval_eta);
            walls_sl(*c,cc->Fifsf,gcval_fifsf);

            comms_off guard(pgc);

            // sigma grid as fnpf_sigma::sigma_ini (the bed is prolonged, not smoothed again)
            FLOOP
            p->sig[FIJK] = p->ZN[KP];

            SLICELOOP4
            {
                k=0;
                p->sig[FIJKm1] = p->ZN[KM1];
                p->sig[FIJKm2] = p->ZN[KM2];
                p->sig[FIJKm3] = p->ZN[KM3];
                k=p->knoz;
                p->sig[FIJKp1] = p->ZN[KP1];
                p->sig[FIJKp2] = p->ZN[KP2];
                p->sig[FIJKp3] = p->ZN[KP3];
            }

            c->pf->fsfdisc_ini(p,cc,pgc,cc->eta,cc->Fifsf);
            c->pf->fsfdisc(p,cc,pgc,cc->eta,cc->Fifsf);
            c->psig->sigma_update(p,cc,pgc,c->pf,cc->eta);
        }

        // Fi
        for(int id : lev[l])
        if(P[id]->fresh)
        prolong_interior_col(*FP(id),sfi);
        fill_col(l,7540+l,sfi);

        for(int id : lev[l])
        if(P[id]->fresh)
        {
            fnpf_amr_patch *c = FP(id);
            comms_off guard(pgc);
            walls_fi(*c,c->c->Fi);
            c->pfu->fsfbc_sig(c->pp,c->c,pgc,c->c->Fifsf,c->c->Fi);
            c->pbu->bedbc_sig(c->pp,c->c,pgc,c->c->Fi,c->pf);
        }
    }
}

// level 0 consistent with the patches
void fnpf_amr::regrid_finish(ghostcell *pgc, int old_total)
{
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;

    restrict_sl(2,[&](int g, int m) -> slice& { return sval(g,2,m); });
    restrict_col([&](int g) -> double* { return gfd(g)->Fi; });

    if(patches_total>0 || old_total>0)
    {
        pgc->gcsl_start4(p0,c0->eta,gcval_eta);
        pgc->gcsl_start4(p0,c0->Fifsf,gcval_fifsf);
        pgc->start7V(p0,c0->Fi,c0->bc,250);
    }
}

// --------------------------------------------------------------------- time stepping
// the tendencies of the patches (kinematic and dynamic FSBC), as fnpf_RK3: before the body
// loads, which need them for the psi_0 data
void fnpf_amr::stage_tendency(lexer *p, fdm_fnpf *c, ghostcell *pgc, int s)
{
    if(!active())
    return;

    double t0 = MPI_Wtime();

    comms_off guard(pgc);

    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        lexer *pp = pc->pp;
        fdm_fnpf *cc = pc->c;

        pp->dt = p->dt;
        pp->dt_old = p->dt_old;
        pp->simtime = p->simtime;
        pp->count = p->count;

        slice4 &ek = *pc->ek, &fk = *pc->fk;
        slice4 &erk1 = *pc->erk1, &erk2 = *pc->erk2;
        slice &Ein = (s==0) ? static_cast<slice&>(cc->eta) : (s==1) ? static_cast<slice&>(erk1) : static_cast<slice&>(erk2);

        lexer *p = pp;
        fdm_fnpf *c = cc;

        pc->pf->kfsfbc(p,c,pgc);
        SLICELOOP4
        ek(i,j) = c->K(i,j);

        pc->pf->dfsfbc(p,c,pgc,Ein);
        SLICELOOP4
        fk(i,j) = c->K(i,j);
    }

    tm[0] += MPI_Wtime()-t0;
}

// the patches form their stage values (as fnpf_RK3) and the footprint of the body, the fine
// values go to the coarser levels, then the cells around the patches are filled
void fnpf_amr::stage_surface(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &Se, slice &Sf, int s)
{
    Se0 = &Se;
    Sf0 = &Sf;
    stg = s;

    if(!active())
    return;

    double t0 = MPI_Wtime();

    for(int id=0; id<(int)P.size(); ++id)
    {
        fnpf_amr_patch *pc = FP(id);

        {
        comms_off guard(pgc);

        fdm_fnpf *cc = pc->c;
        slice4 &ek = *pc->ek, &fk = *pc->fk;
        slice4 &erk1 = *pc->erk1, &erk2 = *pc->erk2, &frk1 = *pc->frk1, &frk2 = *pc->frk2;

        lexer *p = pc->pp;
        fdm_fnpf *c = cc;

        if(s==0)
        {
            SLICELOOP4
            erk1(i,j) = c->eta(i,j) + p->dt*ek(i,j);
            SLICELOOP4
            frk1(i,j) = c->Fifsf(i,j) + p->dt*fk(i,j);
        }
        if(s==1)
        {
            SLICELOOP4
            erk2(i,j) = 0.75*c->eta(i,j) + 0.25*erk1(i,j) + 0.25*p->dt*ek(i,j);
            SLICELOOP4
            frk2(i,j) = 0.75*c->Fifsf(i,j) + 0.25*frk1(i,j) + 0.25*p->dt*fk(i,j);
        }
        if(s==2)
        {
            SLICELOOP4
            c->eta(i,j) = (1.0/3.0)*c->eta(i,j) + (2.0/3.0)*erk2(i,j) + (2.0/3.0)*p->dt*ek(i,j);
            SLICELOOP4
            c->Fifsf(i,j) = (1.0/3.0)*c->Fifsf(i,j) + (2.0/3.0)*frk2(i,j) + (2.0/3.0)*p->dt*fk(i,j);
        }
        }

        if(body!=nullptr)
        body->amr_surface(pc->pp,pgc,id,sval(id,s,0),sval(id,s,1));
    }

    auto sel = [&](int g, int m) -> slice& { return sval(g,s,m); };

    restrict_sl(2,sel);

    pgc->gcsl_start4(p,Se,gcval_eta);
    pgc->gcsl_start4(p,Sf,gcval_fifsf);

    for(int l=1; l<=maxlev; ++l)
    {
        fill_sl(l,2,7500+l,sel);

        for(int id : lev[l])
        {
            walls_sl(*FP(id),sel(id,0),gcval_eta);
            walls_sl(*FP(id),sel(id,1),gcval_fifsf);
        }
    }

    tm[0] += MPI_Wtime()-t0;
}

void fnpf_amr::step_end(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(!active())
    return;

    comms_off guard(pgc);
    for(auto q : P)
    {
        fnpf_amr_patch *pc = FP(q);
        pc->pbu->bedbc_sig(pc->pp,pc->c,pgc,pc->c->Fi,pc->pf);
        pc->pfu->velcalc_sig(pc->pp,pc->c,pgc,pc->c->Fi);
    }
}

// time step of the finest grid, same CFL as fnpf_timestep
void fnpf_amr::timestep(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(maxlev<1)
    return;

    if(p->N48==0)
    return;

    double r=1.0;
    for(auto q : P)
    r = MIN(r, q->pp->DXM/p->DXM);
    r = pgc->globalmin(r);

    if(p->count==0)
    {
        p->dt *= r;
        p->dt_old = p->dt;
        return;
    }

    double depthmax=0.0;
    SLICELOOP4
    depthmax = MAX(depthmax,c->depth(i,j));
    depthmax = pgc->globalmax(depthmax);

    double cu = 1.0e10;
    for(auto q : P)
    {
        lexer *pp = q->pp;
        for(int ii=EXT; ii<EXT+q->nx; ++ii)
        for(int jj=EXT; jj<EXT+q->ny; ++jj)
        if(pp->flagslice4[lij(pp,ii,jj)]>0)
        {
            double dx = MIN(pp->DXN[ii+marge],pp->DYN[jj+marge]);
            cu = MIN(cu, 1.0/((fabs(MAX(p->umax, sqrt(9.81*depthmax)))/dx)));
            cu = MIN(cu, 1.0/((fabs(MAX(p->vmax, sqrt(9.81*depthmax)))/dx)));
        }
    }

    double dtp = pgc->globalmin(p->N47*cu);
    p->dt = MIN(p->dt,dtp);
}
