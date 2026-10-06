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


#include"seastate_amr.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"seastate_exchange.h"
#include"seastate_source.h"
#include"seastate_bathy.h"
#include"seastate_forcing.h"
#include"fdm_seastate.h"
#include"slice4.h"
#include"sliceint4.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<cstdint>
#include<iostream>
#include<iomanip>
#include<mpi.h>
#include<sys/stat.h>
#include<sys/types.h>

namespace
{
// floor(a/2^l) for negative a as well
inline int fsh(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }

inline double mmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;
    return (std::fabs(a)<std::fabs(b)) ? a : b;
}
}

seastate_amr::seastate_amr(lexer *p, ghostcell *pgc) : reefamr(p,pgc), p0(p), nbin(0), nsig(0), fills(0), sweeps(0)
{
    L0 = level0{nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,0.0};
    tm[0]=tm[1]=tm[2]=tm[3]=0.0;
}

seastate_amr::~seastate_amr()
{
    free_patches();
    delete cov0;
}

fdm_seastate* seastate_amr::gfd(int g)
{
    return (g<0) ? L0.e : SP(g)->e;
}

int seastate_amr::level0_covered(int ii, int jj) const
{
    return (cov0!=nullptr) ? (*cov0)(ii,jj) : 0;
}

sliceint4* seastate_amr::gcov(int g)
{
    return (g<0) ? cov0 : SP(g)->cov;
}

// --------------------------------------------------------------------- set-up
void seastate_amr::ini(lexer *p, ghostcell *pgc, const level0 &l0)
{
    L0 = l0;
    nbin = L0.e->grid->nbin;
    nsig = L0.e->grid->nsig;
    zero.assign(nbin,0.0f);

    reefamr_param q;
    q.name = "SEASTATE AMR";
    q.maxlev = p->G1;
    q.regrid = 0;
    q.nbuf = p->G3;
    q.tile = p->G4;
    q.nest = 3;
    q.keep = 0;
    q.ext = 1;
    q.ioband = 0;

    for(int n=0; n<p->G10; ++n)
    {
    q.rbox.push_back(p->G10_xs[n]); q.rbox.push_back(p->G10_xe[n]);
    q.rbox.push_back(p->G10_ys[n]); q.rbox.push_back(p->G10_ye[n]);
    }

    for(int n=0; n<p->G11; ++n)
    {
    q.fbox.push_back(p->G11_xs[n]); q.fbox.push_back(p->G11_xe[n]);
    q.fbox.push_back(p->G11_ys[n]); q.fbox.push_back(p->G11_ye[n]);
    }

    // no refinement next to the sides with boundary spectra or zero gradient (forcing on level 0)
    const double band = double(p->A794)*p->DXM;
    const double xa = p->global_xmin, xb = p->global_xmax, ya = p->global_ymin, yb = p->global_ymax;
    const double big = 1.0e-6*(xb-xa+yb-ya);
    const int sd[4] = {p->A712_xm,p->A712_xp,p->A712_ym,p->A712_yp};

    if(band>0.0)
    for(int s=0; s<4; ++s)
    if(sd[s]!=0)
    {
    double b[4] = {xa-big,xb+big,ya-big,yb+big};
    if(s==0) b[1] = xa+band;
    if(s==1) b[0] = xb-band;
    if(s==2) b[3] = ya+band;
    if(s==3) b[2] = yb-band;
    for(int k=0; k<4; ++k)
    q.fbox.push_back(b[k]);
    }

    configure(q);

    if(maxlev<1)
    return;

    setup(p,pgc);

    // level-0 bed of the whole domain (static), for the interpolated bed of the patches
    {
    vector<double> mine;
    SLICELOOP4
    {
    mine.push_back(double(i+p->origin_i));
    mine.push_back(double(j+p->origin_j));
    mine.push_back(L0.e->bed(i,j));
    }

    const int np = p->mpi_size;
    int nm = int(mine.size());
    vector<int> cnt(np), off(np);
    MPI_Allgather(&nm,1,MPI_INT,cnt.data(),1,MPI_INT,MPI_COMM_WORLD);
    int tot=0;
    for(int r=0; r<np; ++r) {off[r]=tot; tot+=cnt[r];}
    vector<double> all(std::max(tot,1));
    MPI_Allgatherv(mine.data(),nm,MPI_DOUBLE,all.data(),cnt.data(),off.data(),MPI_DOUBLE,MPI_COMM_WORLD);

    bed0.assign(size_t(GNX)*GNY,0.0);
    for(int k=0; k+2<tot; k+=3)
    bed0[size_t(all[k])*GNY + size_t(all[k+1])] = all[k+2];
    }

    cov0 = new sliceint4(p);

    mkdir("./REEF3D_SEASTATE_AMR",0777);

    // initial hierarchy: every pass can add one level
    for(int it=0; it<maxlev; ++it)
    regrid(p,pgc,true);

    if(p->mpirank==0)
    {
    cout<<"SEASTATE AMR: "<<maxlev<<" level(s), "<<patches_total<<" patch(es), "<<cells_total<<" refined cells";
    for(int l=1; l<=maxlev; ++l)
    cout<<(l==1 ? " (" : ", ")<<"level "<<l<<": "<<nlevg[l];
    cout<<(maxlev>0 ? ")" : "")<<endl;

    logout.open("./REEF3D_SEASTATE_AMR/REEF3D_SEASTATE_AMR_log.dat");
    logout<<"# count \t simtime \t patches \t refined cells \t active leaf cells \t iterations \t time sweeps level 0 [s] \t time patches [s] \t time fills [s] \t time restriction [s]"<<endl;
    }

    // memory
    double mb=0.0, act=0.0;
    for(auto c : P)
    {
    mb += double(SP(c)->e->N->bytes() + SP(c)->e->kw->bytes() + SP(c)->e->cg->bytes());
    if(SP(c)->N0!=nullptr)
    mb += double(SP(c)->N0->bytes());
    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    for(int jj=EXT; jj<EXT+c->ny; ++jj)
    if(SP(c)->e->wet(ii,jj)==1 && (*SP(c)->cov)(ii,jj)==0)
    act += 1.0;
    }
    mb = pgc->globalsum(mb)/1048576.0;
    act = pgc->globalsum(act);

    if(p->mpirank==0)
    cout<<"SEASTATE AMR: active leaf cells on the patches "<<long(act)<<", memory of the patches "<<fixed<<setprecision(1)<<mb<<" MB"<<endl<<endl;
    cout.unsetf(ios::floatfield);
    cout<<setprecision(6);
}

reefamr_patch* seastate_amr::patch_new()
{
    return new seastate_amr_patch;
}

// the SEASTATE objects of a patch are built in regrid_static, when the hierarchy is known
void seastate_amr::patch_objects(reefamr_patch*, ghostcell*)
{
}

void seastate_amr::patch_delete(reefamr_patch *q)
{
    seastate_amr_patch *c = SP(q);

    delete c->solv;
    delete c->N0;
    delete c->cov;
    delete c->wU;
    delete c->wD;

    if(c->e!=nullptr)
    {
    delete c->e->N;
    delete c->e->kw;
    delete c->e->cg;
    delete c->e;
    }

    c->e = nullptr;
    c->solv = nullptr;
    c->N0 = nullptr;
    c->cov = nullptr;
    c->wU = c->wD = nullptr;
}

// --------------------------------------------------------------------- bed and environment
// level-l bed: level 0 from the grid, finer levels limited linear from the parent (as SFLOW AMR)
double seastate_amr::bed_at(int l, int I, int J)
{
    if(l==0)
    {
    I = std::max(std::min(I,GNX-1),0);
    J = std::max(std::min(J,GNY-1),0);
    return bed0[size_t(I)*GNY + J];
    }

    const unsigned long long key = (static_cast<unsigned long long>(l)<<58) ^ (static_cast<unsigned long long>(uint32_t(I+(1<<28)))<<29) ^ static_cast<unsigned long long>(uint32_t(J+(1<<28)));
    auto it = bedmemo.find(key);
    if(it!=bedmemo.end())
    return it->second;

    const int Ic = fsh(I,1), Jc = fsh(J,1);
    const double bc = bed_at(l-1,Ic,Jc);
    const int lc = l-1;
    const int nxc = GNX<<lc, nyc = GNY<<lc;

    auto inside = [&](int a, int d) {return a>=0 && a<nxc && d>=0 && d<nyc && flag0(fsh(a,lc),fsh(d,lc))>0;};

    double sx=0.0, sy=0.0;
    if(inside(Ic,Jc) && inside(Ic+1,Jc) && inside(Ic-1,Jc))
    sx = mmod(bed_at(lc,Ic+1,Jc)-bc, bc-bed_at(lc,Ic-1,Jc));
    if(inside(Ic,Jc) && inside(Ic,Jc+1) && inside(Ic,Jc-1))
    sy = mmod(bed_at(lc,Ic,Jc+1)-bc, bc-bed_at(lc,Ic,Jc-1));

    const double ox = (I-2*Ic==0) ? -0.25 : 0.25;
    const double oy = (J-2*Jc==0) ? -0.25 : 0.25;
    const double r = bc + sx*ox + sy*oy;
    bedmemo[key] = r;
    return r;
}

void seastate_amr::environment(seastate_amr_patch &c)
{
    lexer *pp = c.pp;
    fdm_seastate *e = c.e;
    const int l = c.lev;
    const int gnx = GNX<<l, gny = GNY<<l;
    const int m = marge;

    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    {
    const int I = ii-EXT+c.I0, J = jj-EXT+c.J0;
    const bool inside = I>=0 && I<gnx && J>=0 && J<gny;
    const int fl = pp->flagslice4[(ii-pp->imin)*pp->jmax + (jj-pp->jmin)];

    double bed;
    if(L0.bathy!=nullptr)
    bed = L0.bathy->cell(pp->XN[ii+m],pp->XN[ii+1+m],pp->YN[jj+m],pp->YN[jj+1+m],p0->wd-p0->A705);
    else
    bed = bed_at(l,I,J);

    e->bed(ii,jj) = bed;
    e->eta(ii,jj) = 0.0;
    e->U(ii,jj) = 0.0;
    e->V(ii,jj) = 0.0;

        // prescribed current for stand-alone runs (as seastate_f::environment)
        if(p0->A720==1)
        {
        const double xs = p0->A721_xs, xe = p0->A721_xe;
        const double x = pp->XP[ii+m];
        const double w = (xe>xs) ? std::min(std::max((x-xs)/(xe-xs),0.0),1.0) : (x>=xs ? 1.0 : 0.0);
        e->U(ii,jj) = (1.0-w)*p0->A721_us + w*p0->A721_ue;
        }

    e->depth(ii,jj) = (inside && fl>0) ? std::max(p0->wd - bed,0.0) : 0.0;
    e->wet(ii,jj) = (inside && fl>0 && e->depth(ii,jj)>=p0->A705) ? 1 : 0;
    e->wet0(ii,jj) = e->wet(ii,jj);
    }
}

// --------------------------------------------------------------------- refinement criteria
void seastate_amr::tag(int l, vector<unsigned char> &M)
{
    const double dsh = p0->A791;
    const int nco = p0->A792;
    const double rgr = p0->A793;

    if(!(dsh>0.0) && nco<=0 && !(rgr>0.0))
    return;

    auto test = [&](lexer *q, fdm_seastate *e, int ii, int jj, int gi0, int gj0, int gnx, int gny)
    {
        if(e->wet(ii,jj)!=1)
        return false;

        const double d = e->depth(ii,jj);

        if(dsh>0.0 && d<dsh)
        return true;

        if(rgr>0.0)
        {
        const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
            for(int k=0; k<4; ++k)
            {
            const int a = ii+di[k], b = jj+dj[k];
                if(e->wet(a,b)==1)
                {
                const double dn = e->depth(a,b);
                if(std::fabs(dn-d)>rgr*std::max(d,dn))
                return true;
                }
            }
        }

        if(nco>0)
        {
        const int n = std::min(nco,q->margin);
            for(int a=ii-n; a<=ii+n; ++a)
            for(int b=jj-n; b<=jj+n; ++b)
            {
            const int I = a+gi0, J = b+gj0;
            if(I<0 || I>=gnx || J<0 || J>=gny)
            continue;
            if(e->wet(a,b)!=1)
            return true;
            }
        }

        return false;
    };

    if(l==0)
    {
        lexer *p = p0;
        SLICELOOP4
        if(test(p,L0.e,i,j,p->origin_i,p->origin_j,p->gknox,p->gknoy))
        tag_cell(1,i+p->origin_i,j+p->origin_j,nbuf,M);
        return;
    }

    for(int id : lev[l])
    {
    seastate_amr_patch *c = SP(id);
    const int gi0 = c->I0-EXT, gj0 = c->J0-EXT;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(test(c->pp,c->e,ii,jj,gi0,gj0,GNX<<l,GNY<<l))
        tag_cell(l+1,ii+gi0,jj+gj0,nbuf,M);
    }
}

// --------------------------------------------------------------------- regridding
// static data of the new patches: bed, depth, active cells, storage, k, cg, gradients, solver
void seastate_amr::regrid_static(ghostcell*)
{
    for(int id=0; id<(int)P.size(); ++id)
    {
    seastate_amr_patch *c = SP(id);

        if(!c->fresh)
        continue;

    lexer *pp = c->pp;
    c->e = new fdm_seastate(pp);
    c->e->grid = L0.e->grid;

    environment(*c);

    // spectral storage for the active cells of the interior and the ring around it (the inflow of the
    // sweeps); the outer rings of the patch arrays are never used by SEASTATE. Tiles of 2 x 2 cells, so
    // that small patches do not allocate whole 16 x 16 tiles
    vector<int> mask(c->e->wet.V,c->e->wet.V+size_t(pp->imax)*pp->jmax);
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    if(ii<EXT-1 || ii>EXT+c->nx || jj<EXT-1 || jj>EXT+c->ny)
    mask[size_t(ii-pp->imin)*pp->jmax + (jj-pp->jmin)] = 0;

    const int ptile = 2;
    c->e->N  = new seastate_store(pp->imin,pp->jmin,pp->imax,pp->jmax,nbin,ptile);
    c->e->N->build(mask.data());
    c->e->N->fill(0.0f);
    c->e->kw = new seastate_store(pp->imin,pp->jmin,pp->imax,pp->jmax,nsig,ptile);
    c->e->kw->build(mask.data());
    c->e->cg = new seastate_store(pp->imin,pp->jmin,pp->imax,pp->jmax,nsig,ptile);
    c->e->cg->build(mask.data());

    seastate_kinematics(pp,c->e,0.0,c->I0-EXT,c->J0-EXT,GNX<<c->lev,GNY<<c->lev);

        if(p0->A700==1)
        {
        c->N0 = new seastate_store(pp->imin,pp->jmin,pp->imax,pp->jmax,nbin,ptile);
        c->N0->build(mask.data());
        }

    c->cov = new sliceint4(pp);

    c->solv = new seastate_implicit(pp,c->e);
    c->solv->sources(L0.src);
    c->solv->range(EXT,EXT+c->nx-1,EXT,EXT+c->ny-1);

        if(L0.wser!=nullptr)
        {
        c->wU = new slice4(pp);
        c->wD = new slice4(pp);
        c->solv->wind_field(c->wU,c->wD);
        }
    }
}

// fill entries next to the interior, initial spectra of the new patches (from the parent)
void seastate_amr::regrid_state(ghostcell*, vector<reefamr_patch*>&)
{
    for(auto q : P)
    {
    seastate_amr_patch *c = SP(q);
    c->ring.clear();

    c->ringr.clear();

        for(int k=0; k<(int)c->fill.size(); ++k)
        {
        const reefamr_fill &f = c->fill[k];
        const bool inx = f.di>=EXT && f.di<EXT+c->nx;
        const bool iny = f.dj>=EXT && f.dj<EXT+c->ny;

            if((inx && (f.dj==EXT-1 || f.dj==EXT+c->ny)) || (iny && (f.di==EXT-1 || f.di==EXT+c->nx)))
            {
            if(f.kind==2)
            c->ringr.push_back(k);
            else
            c->ring.push_back(k);
            }
        }
    }

    for(int l=1; l<=maxlev; ++l)
    {
        block_down_if(l,nbin,7300+l,
            [&](reefamr_patch *q) {return q->fresh;},
            [&](const reefamr_block &B, int, double *v)
            {
                fdm_seastate *e = gfd(B.g);
                const float *s = (e->wet(B.ic,B.jc)==1) ? e->N->spec(B.ic,B.jc) : nullptr;
                for(int b=0; b<nbin; ++b)
                v[b] = (s!=nullptr) ? double(s[b]) : 0.0;
            },
            [&](reefamr_patch *q, int, int k, const double *v)
            {
                seastate_amr_patch *c = SP(q);
                const int bi = k/(c->ny/2), bj = k%(c->ny/2);
                for(int a=0; a<2; ++a)
                for(int d=0; d<2; ++d)
                {
                    float *s = c->e->N->spec(EXT+2*bi+a,EXT+2*bj+d);
                    if(s!=nullptr && c->e->wet(EXT+2*bi+a,EXT+2*bj+d)==1)
                    for(int b=0; b<nbin; ++b)
                    s[b] = float(v[b]);
                }
            });

        fill_level(l);
    }
}

void seastate_amr::regrid_finish(ghostcell *pgc, int)
{
    covered();

    // restriction without the block plans when every parent cell is on the rank of its patch
    int rem = 0;
    for(auto q : P)
    for(int g : q->rgrid)
    if(g==-3)
    rem = 1;
    remote_blocks = pgc->globalimax(rem)>0;

    // downstream order of the patches of every level for the four quadrants (q 0: +x +y, 1: -x +y,
    // 2: -x -y, 3: +x -y), by the projection of the box centre on the quadrant diagonal
    for(int q=0; q<4; ++q)
    {
    order[q].assign(maxlev+1,vector<int>());
    const double sx = (q==1 || q==2) ? -1.0 : 1.0;
    const double sy = (q==2 || q==3) ? -1.0 : 1.0;

        for(int l=1; l<=maxlev; ++l)
        {
        vector<int> &o = order[q][l];
        o = lev[l];
        std::stable_sort(o.begin(),o.end(),[&](int a, int b)
            {
            const double pa = sx*double(P[a]->I0+P[a]->I1) + sy*double(P[a]->J0+P[a]->J1);
            const double pb = sx*double(P[b]->I0+P[b]->I1) + sy*double(P[b]->J0+P[b]->J1);
            return pa<pb;
            });
        }
    }

    L0.solv->skip(cov0);
    L0.solv->faces(gfaces(-1));

    for(int id=0; id<(int)P.size(); ++id)
    {
    SP(id)->solv->skip(SP(id)->cov);
    SP(id)->solv->faces(gfaces(id));
    }

    restrict_all();
    L0.pex->start(p0,pgc,*L0.e->N);
}

// covered cells of every grid and the faces of the covered cells next to cells solved on their grid
void seastate_amr::covered()
{
    lexer *p = p0;

    IMALOOP
    JMALOOP
    (*cov0)(i,j) = 0;

    for(auto q : P)
    {
    seastate_amr_patch *c = SP(q);
    lexer *pp = c->pp;
    for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
    for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
    (*c->cov)(ii,jj) = 0;
    }

    for(auto q : P)
    {
    seastate_amr_patch *c = SP(q);
    for(size_t k=0; k<c->rgrid.size(); ++k)
    if(c->rgrid[k]>=-1)
    (*gcov(c->rgrid[k]))(c->ric[k],c->rjc[k]) = 1;
    }

    faces.assign(P.size()+1,seastate_amr_faces());
    flinks.assign(maxlev+1,vector<seastate_amr_flink>());

    // the cells a grid solves: level 0 its interior, a patch its box
    auto solved = [&](int g, int a, int d)
    {
        if(g<0)
        return a>=0 && a<p0->knox && d>=0 && d<p0->knoy;
        reefamr_patch *c = P[g];
        return a>=EXT && a<EXT+c->nx && d>=EXT && d<EXT+c->ny;
    };

    const int di[4]={-1,1,0,0}, dj[4]={0,0,-1,1};

    for(int id=0; id<(int)P.size(); ++id)
    {
    seastate_amr_patch *c = SP(id);
    const int nby = c->ny/2;

        for(size_t k=0; k<c->rgrid.size(); ++k)
        {
        const int g = c->rgrid[k];
        if(g<-1)
        continue;

        const int ic = c->ric[k], jc = c->rjc[k];
        fdm_seastate *ec = gfd(g);
        sliceint4 &cv = *gcov(g);
        const int bi = int(k)/nby, bj = int(k)%nby;

            for(int s=0; s<4; ++s)
            {
            const int a = ic+di[s], d = jc+dj[s];

            if(!solved(g,a,d) || ec->wet(a,d)!=1 || cv(a,d)==1)
            continue;

            seastate_amr_flink L;
            L.g = g;
            L.key = seastate_amr_faces::key(ic,jc,s);
            L.id = id;

            const int fi = EXT+2*bi, fj = EXT+2*bj;
            if(s==0) {L.i0=fi;   L.j0=fj;   L.i1=fi;   L.j1=fj+1;}
            if(s==1) {L.i0=fi+1; L.j0=fj;   L.i1=fi+1; L.j1=fj+1;}
            if(s==2) {L.i0=fi;   L.j0=fj;   L.i1=fi+1; L.j1=fj;}
            if(s==3) {L.i0=fi;   L.j0=fj+1; L.i1=fi+1; L.j1=fj+1;}

            gfaces(g)->F[L.key].assign(nbin,0.0f);
            flinks[c->lev].push_back(L);
            }
        }
    }
}

// --------------------------------------------------------------------- composite solve
// the cells next to the interior of one patch from sources on this rank: copy (same level) or constant
// prolongation (coarser level), float to float
void seastate_amr::fill_local(seastate_amr_patch *c)
{
        for(int k : c->ring)
        {
        const reefamr_fill &f = c->fill[k];
        float *d = c->e->N->spec(f.di,f.dj);

        if(d==nullptr)
        continue;

        const float *s = nullptr;
            if(f.kind==0 || f.kind==1)
            {
            fdm_seastate *e = gfd(f.g);
            if(e->wet(f.si,f.sj)==1)
            s = e->N->spec(f.si,f.sj);
            }

        if(s!=nullptr)
        std::copy(s,s+nbin,d);
        else
        std::fill(d,d+nbin,0.0f);
        }
}

void seastate_amr::fill_level(int l)
{
    for(int id : lev[l])
    fill_local(SP(id));

    // cells held by other ranks: the exchange plans of the level
    const reefamr_xplan &X = gplan[l];
    if(X.speer.empty() && X.rpeer.empty())
    {
    ++fills;
    return;
    }

    fill_run_sub(l,nbin,7400+l,
        [&](const reefamr_fill &f, double *v)
        {
            const float *s = nullptr;
            if(f.kind==0 || f.kind==1)
            {
                fdm_seastate *e = gfd(f.g);
                if(e->wet(f.si,f.sj)==1)
                s = e->N->spec(f.si,f.sj);
            }
            for(int b=0; b<nbin; ++b)
            v[b] = (s!=nullptr) ? double(s[b]) : 0.0;
        },
        [&](reefamr_patch *q, int, const reefamr_fill &f, const double *v)
        {
            float *s = SP(q)->e->N->spec(f.di,f.dj);
            if(s!=nullptr)
            for(int b=0; b<nbin; ++b)
            s[b] = float(v[b]);
        },
        [&](reefamr_patch *q) -> const vector<int>* {return &SP(q)->ringr;});

    ++fills;
}

// finest level first: covered cells take the mean of their four children, the faces of covered
// cells next to solved cells the mean of the two fine cells next to the face
void seastate_amr::restrict_all()
{
    for(int l=maxlev; l>=1; --l)
    restrict_level(l);
}

// the level-l patches into their parents
void seastate_amr::restrict_level(int l)
{
    if(remote_blocks)
    {
        block_up(l,nbin,7500+l,
            [&](reefamr_patch *q, int, int k, double *v)
            {
                seastate_amr_patch *c = SP(q);
                const int bi = k/(c->ny/2), bj = k%(c->ny/2);
                for(int b=0; b<nbin; ++b)
                v[b] = 0.0;
                for(int a=0; a<2; ++a)
                for(int d=0; d<2; ++d)
                {
                    const int ii = EXT+2*bi+a, jj = EXT+2*bj+d;
                    const float *s = (c->e->wet(ii,jj)==1) ? c->e->N->spec(ii,jj) : nullptr;
                    if(s!=nullptr)
                    for(int b=0; b<nbin; ++b)
                    v[b] += 0.25*double(s[b]);
                }
            },
            [&](const reefamr_block &B, int, const double *v)
            {
                fdm_seastate *e = gfd(B.g);
                float *s = e->N->spec(B.ic,B.jc);
                if(s!=nullptr && e->wet(B.ic,B.jc)==1)
                for(int b=0; b<nbin; ++b)
                s[b] = float(v[b]);
            });
    }
    else
    for(int id : lev[l])
    {
    seastate_amr_patch *c = SP(id);
    const int nby = c->ny/2;

        for(size_t k=0; k<c->rgrid.size(); ++k)
        {
        const int g = c->rgrid[k];
        if(g<-1)
        continue;

        fdm_seastate *e = gfd(g);
        if(e->wet(c->ric[k],c->rjc[k])!=1)
        continue;

        float *s = e->N->spec(c->ric[k],c->rjc[k]);
        if(s==nullptr)
        continue;

        const int bi = int(k)/nby, bj = int(k)%nby;
        const float *ch[4];
        int n=0;

            for(int a=0; a<2; ++a)
            for(int d=0; d<2; ++d)
            {
            const int ii = EXT+2*bi+a, jj = EXT+2*bj+d;
            if(c->e->wet(ii,jj)==1 && c->e->N->spec(ii,jj)!=nullptr)
            ch[n++] = c->e->N->spec(ii,jj);
            }

            for(int b=0; b<nbin; ++b)
            {
            float v = 0.0f;
            for(int m=0; m<n; ++m)
            v += ch[m][b];
            s[b] = 0.25f*v;
            }
        }
    }

    for(const seastate_amr_flink &L : flinks[l])
    {
    seastate_amr_patch *c = SP(L.id);
    float *F = gfaces(L.g)->F[L.key].data();
    const float *a = (c->e->wet(L.i0,L.j0)==1) ? c->e->N->spec(L.i0,L.j0) : nullptr;
    const float *b = (c->e->wet(L.i1,L.j1)==1) ? c->e->N->spec(L.i1,L.j1) : nullptr;

        for(int n=0; n<nbin; ++n)
        F[n] = 0.5f*((a!=nullptr ? a[n] : 0.0f) + (b!=nullptr ? b[n] : 0.0f));
    }
}

void seastate_amr::step_begin()
{
    for(auto q : P)
    if(SP(q)->N0!=nullptr)
    SP(q)->N0->copy_from(*SP(q)->e->N);
}

void seastate_amr::sweep_level(int l, int q, double rdt, bool refraction, bool fshift)
{
    static const int side0[4] = {0,0,0,0};
    static const vector<float> none;

    for(int id : order[q][l])
    {
    seastate_amr_patch *c = SP(id);
    fill_local(c);
    c->solv->sweep(c->pp,c->e,q,c->N0,rdt,none,side0,refraction,fshift);
    }
}

void seastate_amr::iterate(lexer *p, ghostcell *pgc, const seastate_store *N0, double rdt, const vector<float> &Nb,
                           const int side[4], bool refraction, bool fshift)
{
    for(int q=0; q<4; ++q)
    {
        // down: level 0, then the patches level by level
        double t0 = pgc->timer();
        L0.solv->sweep(p,L0.e,q,N0,rdt,Nb,side,refraction,fshift);
        double t1 = pgc->timer();
        tm[0] += t1-t0;

        for(int l=1; l<=maxlev; ++l)
        {
            t0 = pgc->timer();
            fill_level(l);
            t1 = pgc->timer();
            tm[2] += t1-t0;

            sweep_level(l,q,rdt,refraction,fshift);
            tm[1] += pgc->timer()-t1;
        }

        // up: restriction of the finest level, then the coarser levels again with the outflow of the
        // finer ones, each restricted into its parent
        t0 = pgc->timer();
        restrict_level(maxlev);
        tm[3] += pgc->timer()-t0;

        for(int l=maxlev-1; l>=1; --l)
        {
            t0 = pgc->timer();
            fill_level(l);
            t1 = pgc->timer();
            tm[2] += t1-t0;

            sweep_level(l,q,rdt,refraction,fshift);
            t0 = pgc->timer();
            tm[1] += t0-t1;

            restrict_level(l);
            tm[3] += pgc->timer()-t0;
        }

        // level 0 again, downstream of the level-1 patches of this rank only (upstream cells cannot
        // take anything from the patches within this sweep)
        t0 = pgc->timer();
        if(!lev[1].empty())
        {
        int ilo=p->knox, ihi=-1, jlo=p->knoy, jhi=-1;
            for(int id : lev[1])
            {
            ilo = std::min(ilo,(P[id]->I0>>1)-O0i);
            ihi = std::max(ihi,(P[id]->I1>>1)-O0i);
            jlo = std::min(jlo,(P[id]->J0>>1)-O0j);
            jhi = std::max(jhi,(P[id]->J1>>1)-O0j);
            }
        const bool idown = (q==1 || q==2), jdown = (q==2 || q==3);
        const int ia = idown ? 0 : std::max(ilo-1,0), ib = idown ? std::min(ihi+1,p->knox-1) : p->knox-1;
        const int ja = jdown ? 0 : std::max(jlo-1,0), jb = jdown ? std::min(jhi+1,p->knoy-1) : p->knoy-1;
        L0.solv->range(ia,ib,ja,jb);
        L0.solv->sweep(p,L0.e,q,N0,rdt,Nb,side,refraction,fshift);
        L0.solv->unrange();
        }
        t1 = pgc->timer();
        tm[0] += t1-t0;

        L0.pex->start(p,pgc,*L0.e->N);
        tm[3] += pgc->timer()-t1;

        ++sweeps;
    }
}

// --------------------------------------------------------------------- parameters, wind, output
void seastate_amr::parameters(vector<double> *hs, double &hmax, double &vmin, vector<double> *w)
{
    seastate_param sp;
    const seastate_grid &g = *L0.e->grid;

    for(auto q : P)
    {
    seastate_amr_patch *c = SP(q);
    fdm_seastate *e = c->e;
    lexer *pp = c->pp;

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
        e->Hs(ii,jj)=e->Tm01(ii,jj)=e->Tm10(ii,jj)=e->Tp(ii,jj)=e->dir(ii,jj)=e->spread(ii,jj)=0.0;

        if(e->wet(ii,jj)!=1)
        continue;

        const float *s = e->N->spec(ii,jj);
        sp.compute(g,s);

        e->Hs(ii,jj)     = sp.Hs;
        e->Tm01(ii,jj)   = sp.Tm01;
        e->Tm10(ii,jj)   = sp.Tm10;
        e->Tp(ii,jj)     = sp.Tp;
        e->dir(ii,jj)    = sp.dir;
        e->spread(ii,jj) = sp.spread;

            // leaf cells: interior, not covered
            if(ii>=EXT && ii<EXT+c->nx && jj>=EXT && jj<EXT+c->ny && (*c->cov)(ii,jj)==0)
            {
            if(hs!=nullptr)
            hs->push_back(sp.Hs);
            if(w!=nullptr)
            w->push_back(1.0/double(1<<(2*c->lev)));
            hmax = std::max(hmax,sp.Hs);
            for(int b=0; b<nbin; ++b)
            vmin = std::min(vmin,double(s[b]));
            }
        }
    }
}

void seastate_amr::wind(double tf)
{
    if(L0.wser==nullptr)
    return;

    const int m = marge;

    for(auto q : P)
    {
    seastate_amr_patch *c = SP(q);
    lexer *pp = c->pp;

        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
        double u, v;
        L0.wser->at(tf,pp->XP[ii+m],pp->YP[jj+m],u,v);
        (*c->wU)(ii,jj) = std::sqrt(u*u + v*v);
        (*c->wD)(ii,jj) = std::atan2(v,u);
        }
    }
}

bool seastate_amr::locate(double x, double y, lexer *&q, fdm_seastate *&ee, int &ci, int &cj)
{
    const int m = marge;

    for(int l=maxlev; l>=1; --l)
    for(int id : lev[l])
    {
    seastate_amr_patch *c = SP(id);
    lexer *pp = c->pp;

        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        if(pp->XN[ii+m]<=x && x<pp->XN[ii+1+m])
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        if(pp->YN[jj+m]<=y && y<pp->YN[jj+1+m])
        {
        q = pp;
        ee = c->e;
        ci = ii;
        cj = jj;
        return true;
        }
    }

    return false;
}

// stationary runs: the largest relative change of Hs and the percentage of leaf cells within A 708
void seastate_amr::convergence(lexer *p, int iteration, double change, double percent)
{
    if(p->mpirank!=0)
    return;

    if(!convout.is_open())
    {
    convout.open("./REEF3D_SEASTATE_AMR/REEF3D_SEASTATE_AMR_convergence.dat");
    convout<<"# count \t iteration \t max. relative change of Hs \t percentage of the leaf cell area within A 708"<<endl;
    }

    convout<<p->count<<" \t "<<iteration<<" \t "<<setprecision(6)<<change<<" \t "<<percent<<endl;
}

void seastate_amr::log(lexer *p, ghostcell *pgc, int iterations)
{
    double t[4];
    for(int k=0; k<4; ++k)
    t[k] = pgc->globalmax(tm[k]);

    double act=0.0;
    for(auto q : P)
    {
    seastate_amr_patch *c = SP(q);
    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    for(int jj=EXT; jj<EXT+c->ny; ++jj)
    if(c->e->wet(ii,jj)==1 && (*c->cov)(ii,jj)==0)
    act += 1.0;
    }
    act = pgc->globalsum(act);

    if(p->mpirank==0 && logout.is_open())
    logout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<patches_total<<" \t "<<cells_total<<" \t "<<long(act)<<" \t "<<iterations
          <<" \t "<<setprecision(5)<<t[0]<<" \t "<<t[1]<<" \t "<<t[2]<<" \t "<<t[3]<<endl;
}

// ASCII VTR per patch (interior cells, cell data), a .vtm per print with all patches of all ranks
void seastate_amr::print(lexer *p, ghostcell *pgc)
{
    const int num = printcount_amr;
    const int m = marge;

    for(int id=0; id<(int)P.size(); ++id)
    {
    seastate_amr_patch *c = SP(id);
    lexer *pp = c->pp;
    fdm_seastate *e = c->e;

    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SEASTATE_AMR/REEF3D-SEASTATE-AMR-%08i-%05i-%04i.vtr",num,p->mpirank,id);
    ofstream out(name);

    out<<"<?xml version=\"1.0\"?>"<<endl;
    out<<"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    out<<"<RectilinearGrid WholeExtent=\"0 "<<c->nx<<" 0 "<<c->ny<<" 0 0\">"<<endl;
    out<<"<Piece Extent=\"0 "<<c->nx<<" 0 "<<c->ny<<" 0 0\">"<<endl;
    out<<"<FieldData><DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">"<<setprecision(12)<<p->simtime<<"</DataArray>"
       <<"<DataArray type=\"Int32\" Name=\"level\" NumberOfTuples=\"1\" format=\"ascii\">"<<c->lev<<"</DataArray></FieldData>"<<endl;
    out<<"<Coordinates>"<<endl;
    out<<"<DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">";
    for(int ii=EXT; ii<=EXT+c->nx; ++ii)
    out<<" "<<pp->XN[ii+m];
    out<<"</DataArray>"<<endl;
    out<<"<DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">";
    for(int jj=EXT; jj<=EXT+c->ny; ++jj)
    out<<" "<<pp->YN[jj+m];
    out<<"</DataArray>"<<endl;
    out<<"<DataArray type=\"Float64\" Name=\"z\" format=\"ascii\"> 0</DataArray>"<<endl;
    out<<"</Coordinates>"<<endl;
    out<<"<CellData>"<<endl;

    out<<setprecision(7);
    slice *f[6] = {&e->Hs,&e->Tm01,&e->Tp,&e->dir,&e->spread,&e->depth};
    const char *fn[6] = {"Hs","Tm01","Tp","dir","spread","depth"};

        for(int k=0; k<6; ++k)
        {
        out<<"<DataArray type=\"Float32\" Name=\""<<fn[k]<<"\" format=\"ascii\">";
        for(int jj=EXT; jj<EXT+c->ny; ++jj)
        for(int ii=EXT; ii<EXT+c->nx; ++ii)
        out<<" "<<(*f[k])(ii,jj);
        out<<"</DataArray>"<<endl;
        }

    out<<"<DataArray type=\"Int32\" Name=\"wet\" format=\"ascii\">";
    for(int jj=EXT; jj<EXT+c->ny; ++jj)
    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    out<<" "<<e->wet(ii,jj);
    out<<"</DataArray>"<<endl;

    out<<"<DataArray type=\"Int32\" Name=\"covered\" format=\"ascii\">";
    for(int jj=EXT; jj<EXT+c->ny; ++jj)
    for(int ii=EXT; ii<EXT+c->nx; ++ii)
    out<<" "<<(*c->cov)(ii,jj);
    out<<"</DataArray>"<<endl;

    out<<"</CellData>"<<endl;
    out<<"</Piece>"<<endl<<"</RectilinearGrid>"<<endl<<"</VTKFile>"<<endl;
    }

    // multiblock of all patches
    const int np = p->mpi_size;
    int mine = int(P.size());
    vector<int> cnt(np);
    MPI_Gather(&mine,1,MPI_INT,cnt.data(),1,MPI_INT,0,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
    char name[256];
    snprintf(name,sizeof(name),"./REEF3D_SEASTATE_AMR/REEF3D-SEASTATE-AMR-%08i.vtm",num);
    ofstream out(name);
    out<<"<?xml version=\"1.0\"?>"<<endl;
    out<<"<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">"<<endl;
    out<<"<vtkMultiBlockDataSet>"<<endl;
    int b=0;
        for(int r=0; r<np; ++r)
        for(int id=0; id<cnt[r]; ++id)
        {
        char f[128];
        snprintf(f,sizeof(f),"REEF3D-SEASTATE-AMR-%08i-%05i-%04i.vtr",num,r,id);
        out<<"<DataSet index=\""<<b++<<"\" file=\""<<f<<"\"/>"<<endl;
        }
    out<<"</vtkMultiBlockDataSet>"<<endl<<"</VTKFile>"<<endl;
    }

    ++printcount_amr;
}
