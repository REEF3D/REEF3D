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

#include"seastate_implicit.h"
#include"seastate_exchange.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_dispersion.h"
#include"seastate_source.h"
#include"fdm_seastate.h"
#include"sliceint.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

seastate_implicit::seastate_implicit(lexer *p, fdm_seastate *e) : src(nullptr), wU(nullptr), wD(nullptr), second(false)
{
    for(int k=0; k<4; ++k)
    Nbs[k] = Nbs0[k] = nullptr;

    const seastate_grid &g = *e->grid;

    nsig = g.nsig;
    ndir = g.ndir;

    P.assign(g.nbin,0.0);
    D.assign(g.nbin,0.0);

    m0.assign(4,ndir);
    m1.assign(4,-1);

    for(int m=0; m<ndir; ++m)
    {
    m0[g.quad[m]] = std::min(m0[g.quad[m]],m);
    m1[g.quad[m]] = std::max(m1[g.quad[m]],m);
    }

    Ar.assign(nsig,0.0);
    qc.assign(ndir,0.0);
    qs.assign(ndir,0.0);
    cur.assign(ndir,0.0);
    tAp.assign(ndir,0.0);
    tCp.assign(ndir,0.0);
    tAm.assign(ndir,0.0);
    tCm.assign(ndir,0.0);
    csg.assign(size_t(nsig)*ndir,0.0);
    zero.assign(g.nbin,0.0f);
    zerod.assign(g.nbin,0.0);
    for(vector<double> *v : {&tdi,&tla,&tup,&trh,&tsi,&trd,&tcp,&tdp,&tso})
    v->assign(g.nbin,0.0);
    tok.assign(nsig,1);
    rsum.assign(nsig,0.0);
    wlo.assign(nsig,0);
    whi.assign(nsig,-1);
    thr.assign(nsig,0.0);
    la.assign(ndir,0.0);
    di.assign(ndir,0.0);
    up.assign(ndir,0.0);
    rhs.assign(ndir,0.0);
    sol.assign(ndir,0.0);
    cp.assign(ndir,0.0);
    dp.assign(ndir,0.0);
}

void seastate_implicit::iterate(lexer *p, ghostcell *pgc, fdm_seastate *e, seastate_exchange *pex, const seastate_store *N0,
                                double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    for(int q=0; q<4; ++q)
    {
    sweep(p,e,q,N0,rdt,Nb,side,refraction,fshift);

    if(pex!=nullptr)
    pex->start(p,pgc,*e->N);
    }
}

void seastate_implicit::sweep(lexer *p, fdm_seastate *e, int q, const seastate_store *N0, double rdt,
                              const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    if(m1[q]<m0[q])
    return;

    const bool idown = (q==1 || q==2);
    const bool jdown = (q==2 || q==3);

    // the whole grid, or (mesh refinement) the interior of a patch
    const int ia = ranged ? ri0 : 0, ib = ranged ? ri1 : p->knox-1;
    const int ja = ranged ? rj0 : 0, jb = ranged ? rj1 : p->knoy-1;

    for(int ii=ia; ii<=ib; ++ii)
    for(int jj=ja; jj<=jb; ++jj)
    {
    i = idown ? ia+ib-ii : ii;
    j = jdown ? ja+jb-jj : jj;

        if(skp!=nullptr && (*skp)(i,j)==1)
        {
            // composite sweep (mesh refinement): the finer cells under this one, in its place
            if(vis!=nullptr)
            {
            const int si = i, sj = j;
            (*vis)(si,sj);
            i = si;
            j = sj;
            }
        continue;
        }

        if(e->wet(i,j)==1)
        {
            if(second)
            cell_surfbeat(p,e,q,N0,rdt,Nb,side,refraction);
            else
            cell(p,e,q,N0!=nullptr ? N0->spec(i,j) : nullptr,rdt,Nb,side,refraction,fshift);
        }
    }
}

void seastate_implicit::solve(lexer *p, fdm_seastate *e, int q, int ci, int cj, const seastate_store *N0, double rdt,
                              const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    if(m1[q]<m0[q])
    return;

    i = ci;
    j = cj;

    if(e->wet(i,j)==1)
    cell(p,e,q,N0!=nullptr ? N0->spec(i,j) : nullptr,rdt,Nb,side,refraction,fshift);
}

namespace
{
    using neighbour = seastate_neighbour;

    // fs: side of the neighbour that faces the cell (mesh refinement: spectrum on the face of a covered cell)
    void set_neighbour(lexer *p, fdm_seastate *e, int ni, int nj, const vector<float> &Nb, const vector<float> *const Nbs[4], const int side[4], neighbour &nb,
                       sliceint *skp = nullptr, const seastate_faces *fcs = nullptr, int fs = 0)
    {
        nb = neighbour();

        if(e->wet(ni,nj)==1)
        {
        nb.N  = e->N->spec(ni,nj);

            if(fcs!=nullptr && skp!=nullptr && (*skp)(ni,nj)==1)
            {
            const float *f = fcs->face(ni,nj,fs);
            if(f!=nullptr)
            nb.N = f;
            }

        nb.cg = e->cg->spec(ni,nj);
        nb.U  = e->U(ni,nj);
        nb.V  = e->V(ni,nj);
        return;
        }

        const int gi = ni + p->origin_i;
        const int gj = nj + p->origin_j;

        int s=-1;
        if(gi<0)
        s=0;
        else if(gi>=p->gknox)
        s=1;
        else if(gj<0)
        s=2;
        else if(gj>=p->gknoy)
        s=3;

        if(s>=0 && side[s]==1 && !Nb.empty())
        nb.N = Nb.data();

        // surfbeat (x-) and boundary series (A 711 3): spectrum of the boundary cell
        if(s>=0 && side[s]==1 && Nbs[s]!=nullptr)
        {
        const int k = (s<2) ? nj : ni;
            if(k>=0 && k<((s<2) ? p->knoy : p->knox))
            nb.N = Nbs[s]->data() + size_t(k)*size_t(e->N->nbins());
        }

        if(s>=0 && side[s]==2)
        nb.self = true;
    }

    inline double pos(double a) {return a>0.0 ? a : 0.0;}
    inline double neg(double a) {return a<0.0 ? a : 0.0;}

    // second-order correction of the upwind flux through a face with velocity c:
    // upwind cell u, the cell upstream of it uu, the downstream cell dn (active cells, else 0)
    inline double tvd(double c, const float *uu, const float *u, const float *dn, int b)
    {
        if(uu==nullptr || u==nullptr || dn==nullptr)
        return 0.0;

        const double du = double(dn[b]) - double(u[b]);

        if(!(std::fabs(du)>1.0e-30))
        return 0.0;

        const double r = (double(u[b]) - double(uu[b]))/du;
        const double phi = (r + std::fabs(r))/(1.0 + std::fabs(r));

        return 0.5*c*phi*du;
    }

    inline const float *active(fdm_seastate *e, int i, int j)
    {
        return e->wet(i,j)==1 ? e->N->spec(i,j) : nullptr;
    }
}

// ------------------------------------------------------------------------------------------------
// spectral sparsity (A 795): per cell, frequency and quadrant the range of directions that hold
// energy (rg, lo<<8 | hi, 0xFFFF empty); bins outside the range are zero. A sweep solves per
// frequency the window of the bins that can receive energy: the own range (one more direction on
// each side with refraction, the next frequencies with frequency shift, the lower frequencies of the
// triads), the ranges of the four neighbours, the ends of the quadrant next to energy in the
// neighbouring quadrants; the whole quadrant where wind input or DIA can fill any bin.

// total energy of a spectrum
double seastate_implicit::energy(const seastate_grid &g, const float *N) const
{
    double et = 0.0;
    for(int l=0; l<nsig; ++l)
    {
    double s = 0.0;
    const float *Nl = N + g.bin(l,0);
    for(int m=0; m<ndir; ++m)
    s += double(Nl[m]);
    et += s*g.sig[l]*g.dsig[l]*g.dtheta;
    }
    return et;
}

// the cell holds its active ranges: solved by this solver (interior or patch box), wet, not covered
bool seastate_implicit::managed(lexer *p, fdm_seastate *e, int ci, int cj) const
{
    const int ia = ranged ? ri0 : 0, ib = ranged ? ri1 : p->knox-1;
    const int ja = ranged ? rj0 : 0, jb = ranged ? rj1 : p->knoy-1;

    if(ci<ia || ci>ib || cj<ja || cj>jb)
    return false;

    if(e->wet(ci,cj)!=1 || e->N->spec(ci,cj)==nullptr)
    return false;

    return skp==nullptr || (*skp)(ci,cj)==0;
}

uint16_t *seastate_implicit::ranges(lexer *p, fdm_seastate *e, int ci, int cj)
{
    if(rg.empty())
    {
    rni = p->imax; rnj = p->jmax; rimin = p->imin; rjmin = p->jmin;
    rg.assign(size_t(rni)*rnj*4*nsig,0xFFFF);
    rv.assign(size_t(rni)*rnj,uint16_t(0));
    rb.assign(size_t(rni)*rnj*4,uint16_t(0xFFFF));
    }

    const size_t c = size_t(ci-rimin)*rnj + (cj-rjmin);
    uint16_t *r = &rg[c*4*nsig];

    if(rv[c]==0)
    {
    // first visit: ranges of the present spectrum (all quadrants), zero below the threshold
    const seastate_grid &g = *e->grid;
    float *N = e->N->spec(ci,cj);
    const double et = energy(g,N);

        for(int qq=0; qq<4; ++qq)
        {
        const int ma = m0[qq], nq = m1[qq]-m0[qq]+1;
            for(int l=0; l<nsig; ++l)
            {
            float *Nl = N + g.bin(l,ma);
            const double t = eps*et/(g.sig[l]*g.dsig[l]*g.dtheta);
            int lo = nq, hi = -1;
                for(int n=0; n<nq; ++n)
                if(double(Nl[n])>t)
                {
                lo = std::min(lo,n);
                hi = n;
                }
                for(int n=0; n<nq; ++n)
                if(n<lo || n>hi)
                Nl[n] = 0.0f;
            r[qq*nsig+l] = (hi<0) ? 0xFFFF : uint16_t((lo<<8) | hi);
            }
        }

    rv[c] = 1 | 2;
    for(int qq=0; qq<4; ++qq)
    summary(ci,cj,qq,r);
    }

    return r;
}

// quadrant summary of a cell: energy in the quadrant, in its first and its last direction
void seastate_implicit::summary(int ci, int cj, int q, const uint16_t *r)
{
    const size_t c = size_t(ci-rimin)*rnj + (cj-rjmin);
    const int nq = m1[q]-m0[q]+1;
    bool any=false, lo=false, hi=false;
    int la=nsig, lb=-1;

    for(int l=0; l<nsig; ++l)
    {
    const uint16_t rr = r[q*nsig+l];
    if(rr==0xFFFF)
    continue;
    any = true;
    la = std::min(la,l);
    lb = l;
    if((rr>>8)==0) lo = true;
    if(int(rr&0xFF)==nq-1) hi = true;
    }

    uint16_t v = rv[c] & uint16_t(~((4<<q) | (64<<(2*q)) | (128<<(2*q))));
    if(any) v |= uint16_t(4<<q);
    if(lo)  v |= uint16_t(64<<(2*q));
    if(hi)  v |= uint16_t(128<<(2*q));
    rv[c] = v;
    rb[c*4+q] = any ? uint16_t((la<<8) | lb) : uint16_t(0xFFFF);
}

static inline void hull(int &lo, int &hi, uint16_t r)
{
    if(r==0xFFFF)
    return;
    lo = std::min(lo,int(r>>8));
    hi = std::max(hi,int(r&0xFF));
}

void seastate_implicit::windows(lexer *p, fdm_seastate *e, int q, int ci, int cj, const float *N,
                                const neighbour &W, const neighbour &E, const neighbour &S, const neighbour &Nn,
                                bool rf, bool fs)
{
    const seastate_grid &g = *e->grid;
    const int ma = m0[q], nq = m1[q]-m0[q]+1;
    const uint16_t *r = ranges(p,e,ci,cj);

    // the whole quadrant where the source terms can fill any bin
    if(src!=nullptr && src->fills_spectrum())
    {
    for(int l=0; l<nsig; ++l)
    {
    wlo[l] = 0;
    whi[l] = nq-1;
    }
    return;
    }

    const int qm = (q+3)%4, qp = (q+1)%4;
    const int nqm = m1[qm]-m0[qm]+1;
    const uint16_t sc = rv[size_t(ci-rimin)*rnj + (cj-rjmin)];
    const int tr = (src!=nullptr && (sc & 2)) ? src->triad_reach() : 0;

    // neighbours: their ranges if this solver holds them, else the positive bins of their spectrum
    const neighbour *nb[4] = {&W,&E,&S,&Nn};
    const int nbi[4] = {ci-1,ci+1,ci,ci}, nbj[4] = {cj,cj,cj-1,cj+1};
    const uint16_t *nr[4] = {nullptr,nullptr,nullptr,nullptr};

    bool maybe = (sc & (4<<q)) || (rf && ((sc & (128<<(2*qm))) || (sc & (64<<(2*qp)))));

    for(int k=0; k<4; ++k)
    if(nb[k]->N!=nullptr)
    {
        if(nb[k]->cg!=nullptr && managed(p,e,nbi[k],nbj[k]) && nb[k]->N==e->N->spec(nbi[k],nbj[k]))
        {
        nr[k] = ranges(p,e,nbi[k],nbj[k]);
        if(rv[size_t(nbi[k]-rimin)*rnj + (nbj[k]-rjmin)] & (4<<q))
        maybe = true;
        }
        else
        maybe = true;
    }

    for(int l=0; l<nsig; ++l)
    {
    wlo[l] = 0;
    whi[l] = -1;
    }

    // nothing in this quadrant, its neighbours or (refraction) at the ends of the next quadrants
    if(!maybe)
    return;

    // band of frequencies that can receive energy
    const size_t cc = size_t(ci-rimin)*rnj + (cj-rjmin);
    int La = nsig, Lb = -1;
    auto band = [&](uint16_t b, int dn, int up)
    {
        if(b==0xFFFF)
        return;
        La = std::min(La,std::max(int(b>>8)-dn,0));
        Lb = std::max(Lb,std::min(int(b&0xFF)+up,nsig-1));
    };

    band(rb[cc*4+q],fs ? 1 : 0,(fs ? 1 : 0) - tr);
    if(rf)
    {
    if(sc & (128<<(2*qm))) band(rb[cc*4+qm],0,0);
    if(sc & (64<<(2*qp)))  band(rb[cc*4+qp],0,0);
    }
    for(int k=0; k<4; ++k)
    if(nb[k]->N!=nullptr)
    {
        if(nr[k]!=nullptr)
        band(rb[(size_t(nbi[k]-rimin)*rnj + (nbj[k]-rjmin))*4+q],0,0);
        else
        {
        La = 0;
        Lb = nsig-1;
        }
    }

    for(int l=La; l<=Lb; ++l)
    {
    int lo = nq, hi = -1;

    const uint16_t own = r[q*nsig+l];
    hull(lo,hi,own);

        // refraction: one more direction on each side, and the quadrant ends next to energy in the
        // neighbouring quadrants
        if(rf)
        {
            if(own!=0xFFFF)
            {
            lo = std::max(lo-1,0);
            hi = std::min(hi+1,nq-1);
            }

        const uint16_t rm = r[qm*nsig+l], rp = r[qp*nsig+l];
        if(rm!=0xFFFF && int(rm&0xFF)==nqm-1)
        {lo = std::min(lo,0); hi = std::max(hi,0);}
        if(rp!=0xFFFF && int(rp>>8)==0)
        {lo = std::min(lo,nq-1); hi = std::max(hi,nq-1);}
        }

        // frequency shift: the next frequencies
        if(fs)
        {
        if(l>0)      hull(lo,hi,r[q*nsig+l-1]);
        if(l<nsig-1) hull(lo,hi,r[q*nsig+l+1]);
        }

        // triads: the lower frequencies that feed this one
        for(int k=l+tr; k<l; ++k)
        if(k>=0)
        hull(lo,hi,r[q*nsig+k]);

        // inflow from the neighbours
        for(int k=0; k<4; ++k)
        if(nb[k]->N!=nullptr)
        {
            if(nr[k]!=nullptr)
            hull(lo,hi,nr[k][q*nsig+l]);
            else
            {
            const float *Nl = nb[k]->N + g.bin(l,ma);
            int n=0;
            while(n<nq && !(Nl[n]>0.0f)) ++n;
            if(n<nq)
            {
            int m=nq-1;
            while(!(Nl[m]>0.0f)) --m;
            lo = std::min(lo,n);
            hi = std::max(hi,m);
            }
            }
        }

    wlo[l] = lo;
    whi[l] = hi;
    }
}

// after the solve: the new ranges of the solved windows, zero below the threshold at their ends
void seastate_implicit::keep(lexer *p, fdm_seastate *e, int q, int ci, int cj)
{
    const seastate_grid &g = *e->grid;
    const int ma = m0[q];
    uint16_t *r = ranges(p,e,ci,cj);
    float *N = e->N->spec(ci,cj);

    for(int l=lmin; l<=lmax; ++l)
    {
    const int na = wlo[l], nb = whi[l];
    if(na>nb)
    {
    r[q*nsig+l] = 0xFFFF;
    continue;
    }

    float *Nl = N + g.bin(l,ma);
    int lo = nb+1, hi = na-1;

        for(int n=na; n<=nb; ++n)
        if(double(Nl[n])>thr[l])
        {
        if(lo>nb) lo = n;
        hi = n;
        }

        for(int n=na; n<=nb; ++n)
        if(n<lo || n>hi)
        Nl[n] = 0.0f;

    r[q*nsig+l] = (hi<lo) ? 0xFFFF : uint16_t((lo<<8) | hi);
    }

    summary(ci,cj,q,r);
}

// second-order fluxes: spectrum of an active cell within the solved area or the ring around it,
// not covered by a finer grid
const float *seastate_implicit::usable(lexer *p, fdm_seastate *e, int ci, int cj) const
{
    const int ia = ranged ? ri0-1 : -2, ib = ranged ? ri1+1 : p->knox+1;
    const int ja = ranged ? rj0-1 : -2, jb = ranged ? rj1+1 : p->knoy+1;

    if(ci<ia || ci>ib || cj<ja || cj>jb)
    return nullptr;

    if(e->wet(ci,cj)!=1)
    return nullptr;

    if(skp!=nullptr && (*skp)(ci,cj)==1)
    return nullptr;

    return e->N->spec(ci,cj);
}

// solution of frequency l (directions wlo..whi of the quadrant starting at ma) into N: zero after
// a vanishing pivot, the action density limiter of the source terms (|change per iteration| <= dmax,
// dmax < 0: none), clipped at zero with the second-order fluxes; with spectral sparsity the bins below
// the threshold at the ends of the window are set to zero (the active range of the frequency)
void seastate_implicit::store(float *N, int l, int ma, int nq, double dmax, bool clip)
{
    const int na = wlo[l], nb = whi[l], o = l*ndir;
    float *Nl = N + l*ndir + ma;

    if(na>nb)
    return;

    if(tok[l]==0)
    {
    for(int n=na; n<=nb; ++n)
    Nl[n] = 0.0f;
    return;
    }

    if(dmax>=0.0)
    for(int n=na; n<=nb; ++n)
    {
    const double v = double(Nl[n]);
    tso[o+n] = std::min(std::max(tso[o+n],v-dmax),v+dmax);
    }

    if(clip)
    for(int n=na; n<=nb; ++n)
    tso[o+n] = std::max(tso[o+n],0.0);

    for(int n=na; n<=nb; ++n)
    Nl[n] = float(tso[o+n]);
}

void seastate_implicit::cell(lexer *p, fdm_seastate *e, int q, const float *N0, double rdt,
                             const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    const seastate_grid &g = *e->grid;

    float *N = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cgc = e->cg->spec(i,j);

    const int ic = i, jc = j;

    neighbour W,E,S,Nn;
    set_neighbour(p,e,ic-1,jc,Nb,Nbs,side,W,skp,fcs,1);
    set_neighbour(p,e,ic+1,jc,Nb,Nbs,side,E,skp,fcs,0);
    set_neighbour(p,e,ic,jc-1,Nb,Nbs,side,S,skp,fcs,3);
    set_neighbour(p,e,ic,jc+1,Nb,Nbs,side,Nn,skp,fcs,2);

    const double rdx = 1.0/p->DXN[IP];
    const double rdy = 1.0/p->DYN[JP];
    const double rdth = 1.0/g.dtheta;

    const double d = e->depth(ic,jc);
    const double U = e->U(ic,jc), V = e->V(ic,jc);
    const bool kin = (e->refr(ic,jc)==1);

    const double ddx = e->ddx(ic,jc), ddy = e->ddy(ic,jc);
    const double dUdx = e->dUdx(ic,jc), dUdy = e->dUdy(ic,jc), dVdx = e->dVdx(ic,jc), dVdy = e->dVdy(ic,jc);
    const double dddt = e->dddt(ic,jc);

    const int ma = m0[q], mb = m1[q], nq = mb-ma+1;

    // frequency shift only where c_sigma can be nonzero (depth change in time or along the current,
    // current gradients); without it the frequencies are independent
    const double dep = dddt + U*ddx + V*ddy;
    const bool fs = fshift && kin && (dep!=0.0 || dUdx!=0.0 || dUdy!=0.0 || dVdx!=0.0 || dVdy!=0.0);
    const bool rf = refraction && kin;

    // wind field: the wind of the cell
    if(src!=nullptr && wU!=nullptr)
    src->set_wind((*wU)(ic,jc),(*wD)(ic,jc));

    // directions solved per frequency: the whole quadrant, or (spectral sparsity, A 795) the window
    // of the bins that can hold energy
    const bool sparse = (eps>0.0);
    int wmin = 0, wmax = nq-1;

    if(sparse)
    {
    windows(p,e,q,ic,jc,N,W,E,S,Nn,rf,fs);

    wmin = nq; wmax = -1;
    lmin = nsig; lmax = -1;
        for(int l=0; l<nsig; ++l)
        if(wlo[l]<=whi[l])
        {
        wmin = std::min(wmin,wlo[l]);
        wmax = std::max(wmax,whi[l]);
        lmin = std::min(lmin,l);
        lmax = l;
        }

        // nothing can arrive in this quadrant
        if(wmax<0)
        {
        i = ic;
        j = jc;
        return;
        }
    }
    else
    {
    lmin = 0;
    lmax = nsig-1;
        for(int l=0; l<nsig; ++l)
        {
        wlo[l] = 0;
        whi[l] = nq-1;
        }
    }

    // maximum energy (A 737 1, SWAN SINTGRL), then the source terms from the latest spectrum,
    // for the directions that are solved; with spectral sparsity the row sums over the active ranges
    if(src!=nullptr)
    {
        if(sparse)
        {
        const uint16_t *r = ranges(p,e,ic,jc);
            for(int l=0; l<nsig; ++l)
            {
            double s = 0.0;
                for(int qq=0; qq<4; ++qq)
                {
                const uint16_t rr = r[qq*nsig+l];
                if(rr==0xFFFF)
                continue;
                const float *Nl = N + g.bin(l,m0[qq]);
                for(int n=int(rr>>8); n<=int(rr&0xFF); ++n)
                s += double(Nl[n]);
                }
            rsum[l] = s;
            }
        src->set_rows(rsum.data());
        }

    src->cap(N,d);
    src->compute(N,d,kc,cgc,P.data(),D.data(),ma+wmin,ma+wmax,lmin,lmax);
    }

    const double *lim = (src!=nullptr) ? src->limit() : nullptr;

    // triads of this cell (the window of the next sweeps includes the lower frequencies that feed a bin)
    if(sparse && src!=nullptr)
    {
    const size_t c = size_t(ic-rimin)*rnj + (jc-rjmin);
    rv[c] = src->triads_active() ? (rv[c] | 2) : (rv[c] & ~2);
    }

    // threshold of the sparse spectrum: bins with E < A 795 E_tot of the cell are set to zero
    if(sparse)
    {
    const double et = (src!=nullptr) ? src->Etot : energy(g,N);
    for(int l=0; l<nsig; ++l)
    thr[l] = eps*et/(g.sig[l]*g.dsig[l]*g.dtheta);
    }

    // coefficients of the cell that do not depend on the bin: per direction of the quadrant (current
    // shear, refraction factors of the theta faces) and per frequency (depth refraction
    // sig/sinh(2kd), c_sigma); then one pass over the bins

    for(int n=0; n<nq; ++n)
    {
    const int m = ma+n;
    const double cs = g.costh[m], sn = g.sinth[m];

    qc[n] = cs;
    qs[n] = sn;
    cur[n] = cs*cs*dUdx + sn*cs*(dUdy+dVdx) + sn*sn*dVdy;

        if(!rf)
        tAp[n] = tCp[n] = tAm[n] = tCm[n] = 0.0;
        else
        {
        const int mm = (m-1+ndir)%ndir;
        const double sp = g.sinthf[m],  cp_ = g.costhf[m];
        const double sm = g.sinthf[mm], cm  = g.costhf[mm];

        tAp[n] = sp*ddx - cp_*ddy;
        tCp[n] = sp*cp_*(dUdx-dVdy) + sp*sp*dVdx - cp_*cp_*dUdy;
        tAm[n] = sm*ddx - cm*ddy;
        tCm[n] = sm*cm*(dUdx-dVdy) + sm*sm*dVdx - cm*cm*dUdy;
        }
    }

    if(kin)
    for(int l=0; l<nsig; ++l)
    if(fs || wlo[l]<=whi[l])
    Ar[l] = seastate_refraction(g.sig[l],kc[l],d);

    // c_sigma at the cell centre, per frequency and direction of the quadrant (zero without
    // frequency shift)
    if(!fs)
    for(int l=0; l<nsig; ++l)
    for(int n=0; n<nq; ++n)
    csg[size_t(l)*ndir+n] = 0.0;
    else
    for(int l=0; l<nsig; ++l)
    {
    const double a = double(kc[l])*Ar[l]*dep;
    const double bk = double(cgc[l])*double(kc[l]);
    double *c = &csg[size_t(l)*ndir];

        for(int n=0; n<nq; ++n)
        c[n] = a - bk*cur[n];
    }

    // neighbour spectra: zero spectrum where there is no inflow
    const float *NW = W.N ? W.N : zero.data();
    const float *NE = E.N ? E.N : zero.data();
    const float *NS = S.N ? S.N : zero.data();
    const float *NN = Nn.N ? Nn.N : zero.data();
    const bool self = W.self || E.self || S.self || Nn.self;
    const float *Nold = (N0!=nullptr) ? N0 : N;
    const double *Dsrc = (src!=nullptr) ? D.data() : zerod.data();
    const double *Psrc = (src!=nullptr) ? P.data() : zerod.data();

    // second-order geographic fluxes (A 796 2): the cells up to two away (active, solved or the ring
    // around a patch, not covered by a finer grid)
    const bool so = (order2 && second==false);
    const float *aW=nullptr, *aE=nullptr, *aS=nullptr, *aN=nullptr, *aWW=nullptr, *aEE=nullptr, *aSS=nullptr, *aNN=nullptr;

    if(so)
    {
    aW  = usable(p,e,ic-1,jc); aE  = usable(p,e,ic+1,jc); aS  = usable(p,e,ic,jc-1); aN  = usable(p,e,ic,jc+1);
    aWW = usable(p,e,ic-2,jc); aEE = usable(p,e,ic+2,jc); aSS = usable(p,e,ic,jc-2); aNN = usable(p,e,ic,jc+2);
    }

    // A: the tridiagonal systems in theta of all frequencies (arrays [l*ndir + n], directions
    // wlo..whi of each frequency); the inflow from the next lower frequency (c_sigma > 0) is added
    // in C with its latest value
    auto assemble = [&](int l, int na, int nb)
    {
    const double cgl = cgc[l];
    const double A = kin ? Ar[l] : 0.0;
    const double rdsig = 1.0/g.dsig[l];

    // geographic face velocities: (c_cell + c_neighbour)/2 = cg_face cos(theta) + U_face
    const double cgW = W.cg  ? 0.5*(cgl + double(W.cg[l]))  : cgl, UW = W.cg  ? 0.5*(U + W.U)  : U;
    const double cgE = E.cg  ? 0.5*(cgl + double(E.cg[l]))  : cgl, UE = E.cg  ? 0.5*(U + E.U)  : U;
    const double cgS = S.cg  ? 0.5*(cgl + double(S.cg[l]))  : cgl, VS = S.cg  ? 0.5*(V + S.V)  : V;
    const double cgN = Nn.cg ? 0.5*(cgl + double(Nn.cg[l])) : cgl, VN = Nn.cg ? 0.5*(V + Nn.V) : V;

    const double *cc = &csg[size_t(l)*ndir];
    const double *cm = (l>0)      ? &csg[size_t(l-1)*ndir] : cc;
    const double *cu = (l<nsig-1) ? &csg[size_t(l+1)*ndir] : cc;

    const int b0 = g.bin(l,ma);
    const float *Nlp = (l<nsig-1) ? N + g.bin(l+1,ma) : zero.data();
    const int o = l*ndir;

    const double *__restrict Dl = Dsrc + b0, *__restrict Pl = Psrc + b0;
    const float *__restrict Nol = Nold + b0, *__restrict NWl = NW + b0, *__restrict NEl = NE + b0;
    const float *__restrict NSl = NS + b0, *__restrict NNl = NN + b0, *__restrict Nlp_ = Nlp;
    const double *__restrict qc_ = qc.data(), *__restrict qs_ = qs.data();
    const double *__restrict tAp_ = tAp.data(), *__restrict tCp_ = tCp.data(), *__restrict tAm_ = tAm.data(), *__restrict tCm_ = tCm.data();
    const double *__restrict cc_ = cc, *__restrict cu_ = cu, *__restrict cm_ = cm;
    double *__restrict di_ = &tdi[o], *__restrict la_ = &tla[o], *__restrict up_ = &tup[o];
    double *__restrict rh_ = &trh[o], *__restrict si_ = &tsi[o];

        // branch-free, vectorised over the directions
        #pragma GCC ivdep
        for(int n=na; n<=nb; ++n)
        {
        const double cs = qc_[n], sn = qs_[n];

        // c_sigma at the faces l+1/2 and l-1/2 (zero without frequency shift; range ends: outflow only)
        const double csc = cc_[n];
        const double csu = 0.5*(csc + cu_[n]);
        const double csl = 0.5*(csc + cm_[n]);

        // c_theta at the faces m+1/2 and m-1/2 (zero without refraction)
        const double ctp = A*tAp_[n] + tCp_[n];
        const double ctm = A*tAm_[n] + tCm_[n];

        const double cxw = cgW*cs + UW, cxe = cgE*cs + UE;
        const double cys = cgS*sn + VS, cyn = cgN*sn + VN;

        // diagonal: 1/dt + outflow (+ D)
        di_[n] = rdt + (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy
                     + (pos(csu) - neg(csl))*rdsig + Dl[n] + (pos(ctp) - neg(ctm))*rdth;

        // right-hand side: old time level + inflow from the neighbours and from l+1 (+ P)
        rh_[n] = rdt*double(Nol[n]) + Pl[n]
               + pos(cxw)*rdx*double(NWl[n]) - neg(cxe)*rdx*double(NEl[n])
               + pos(cys)*rdy*double(NSl[n]) - neg(cyn)*rdy*double(NNl[n])
               - neg(csu)*rdsig*double(Nlp_[n]);

        // theta: within the window implicit (tridiagonal)
        la_[n] = -pos(ctm)*rdth;
        up_[n] =  neg(ctp)*rdth;
        si_[n] = pos(csl)*rdsig;
        }

    // theta: outside the window (other quadrants, empty bins) with the latest values
    {
    const double ctm0 = A*tAm[na] + tCm[na], ctp1 = A*tAp[nb] + tCp[nb];
    rh_[na] += pos(ctm0)*rdth*double(N[g.bin(l,(ma+na-1+ndir)%ndir)]);
    rh_[nb] -= neg(ctp1)*rdth*double(N[g.bin(l,(ma+nb+1)%ndir)]);
    la_[na] = 0.0;
    up_[nb] = 0.0;
    }

        // zero-gradient sides: inflow of the cell's own spectrum, implicit while the non-theta part
        // of the diagonal stays at least half of its value (M-matrix), otherwise with the latest value
        if(self)
        for(int n=na; n<=nb; ++n)
        {
        const double cs = qc[n], sn = qs[n];
        const double cxw = cgW*cs + UW, cxe = cgE*cs + UE;
        const double cys = cgS*sn + VS, cyn = cgN*sn + VN;
        const double ctp = A*tAp[n] + tCp[n], ctm = A*tAm[n] + tCm[n];
        const double dg = di_[n] - (pos(ctp) - neg(ctm))*rdth;

        double aself = 0.0;
        if(W.self)  aself += pos(cxw)*rdx;
        if(E.self)  aself -= neg(cxe)*rdx;
        if(S.self)  aself += pos(cys)*rdy;
        if(Nn.self) aself -= neg(cyn)*rdy;

            if(aself>0.0)
            {
                if(dg-aself>=0.5*dg)
                di_[n] -= aself;
                else
                rh_[n] += aself*double(N[b0+n]);
            }
        }

        // second-order upwind geographic fluxes (A 796 2, as SWAN SORDUP): the face value is extrapolated
        // from the two upwind cells, F = c (3/2 N_up - 1/2 N_upup). Implicit in the cell, the upwind cells
        // are solved before it in the sweep; faces without two active upwind cells stay first order. Not
        // monotone: the solution is clipped at zero (store)
        if(so)
        for(int n=na; n<=nb; ++n)
        {
        const int b = b0+n;
        const double cs = qc[n], sn = qs[n];
        const double cxw = cgW*cs + UW, cxe = cgE*cs + UE;
        const double cys = cgS*sn + VS, cyn = cgN*sn + VN;
        double dd = 0.0, rr = 0.0;

        // x+ face
        if(cxe>0.0 && aW)        {dd += 0.5*cxe*rdx;  rr += 0.5*cxe*rdx*double(aW[b]);}
        if(cxe<0.0 && aE && aEE) {rr -= 0.5*cxe*rdx*(double(aE[b]) - double(aEE[b]));}
        // x- face
        if(cxw>0.0 && aW && aWW) {rr += 0.5*cxw*rdx*(double(aW[b]) - double(aWW[b]));}
        if(cxw<0.0 && aE)        {dd -= 0.5*cxw*rdx;  rr -= 0.5*cxw*rdx*double(aE[b]);}
        // y+ face
        if(cyn>0.0 && aS)        {dd += 0.5*cyn*rdy;  rr += 0.5*cyn*rdy*double(aS[b]);}
        if(cyn<0.0 && aN && aNN) {rr -= 0.5*cyn*rdy*(double(aN[b]) - double(aNN[b]));}
        // y- face
        if(cys>0.0 && aS && aSS) {rr += 0.5*cys*rdy*(double(aS[b]) - double(aSS[b]));}
        if(cys<0.0 && aN)        {dd -= 0.5*cys*rdy;  rr -= 0.5*cys*rdy*double(aN[b]);}

        di_[n] += dd;
        rh_[n] += rr;
        }
    };

    for(int l=lmin; l<=lmax; ++l)
    if(wlo[l]<=whi[l])
    assemble(l,wlo[l],whi[l]);

    // B: elimination factors of all frequencies (independent, the inner loop over l hides the
    // latency of the divisions); a frequency with a vanishing pivot is set to zero
    // (with spectral sparsity the windows are short: one frequency after the other, in C)
    if(!sparse)
    {
    for(int l=lmin; l<=lmax; ++l)
    tok[l] = 1;

    for(int n=wmin; n<=wmax; ++n)
    for(int l=lmin; l<=lmax; ++l)
    {
    if(n<wlo[l] || n>whi[l])
    continue;

    const int x = l*ndir + n;
    const double den = tdi[x] - (n>wlo[l] ? tla[x]*tcp[x-1] : 0.0);

        if(!(den>1.0e-300))
        {
        tok[l] = 0;
        trd[x] = tcp[x] = 0.0;
        continue;
        }

    trd[x] = 1.0/den;
    tcp[x] = tup[x]*trd[x];
    }
    }

    // one frequency on its own: factors, the inflow from l-1 (frequency shift), substitution
    auto solve1 = [&](int l)
    {
        const int o = l*ndir, na = wlo[l], nb = whi[l];
        if(na>nb)
        return;

        tok[l] = 1;
        for(int n=na; n<=nb; ++n)
        {
        const int x = o+n;
        const double den = tdi[x] - (n>na ? tla[x]*tcp[x-1] : 0.0);
            if(!(den>1.0e-300))
            {
            tok[l] = 0;
            return;
            }
        trd[x] = 1.0/den;
        tcp[x] = tup[x]*trd[x];
        }

        if(fs && l>0)
        {
        const float *Nlm = N + g.bin(l-1,ma);
        for(int n=na; n<=nb; ++n)
        trh[o+n] += tsi[o+n]*double(Nlm[n]);
        }

        for(int n=na; n<=nb; ++n)
        tdp[o+n] = (trh[o+n] - (n>na ? tla[o+n]*tdp[o+n-1] : 0.0))*trd[o+n];

        tso[o+nb] = tdp[o+nb];
        for(int n=nb-1; n>=na; --n)
        tso[o+n] = tdp[o+n] - tcp[o+n]*tso[o+n+1];
    };

    // spectral sparsity: where the solution at an end of the window is above the threshold (energy
    // refracted further than the window within this solve), the window grows and the frequency is
    // solved again, with the source terms of the new bins
    auto widen = [&](int l)
    {
        for(int tries=0; tries<8; ++tries)
        {
        const int na = wlo[l], nb = whi[l];
        if(na>nb || tok[l]==0)
        return;

        // estimate of the bin beyond each end: its inflow through the theta face from the end bin
        // over the diagonal of the end bin
        const float *Nl = N + g.bin(l,ma);
        const double A = kin ? Ar[l] : 0.0;
        const double ctm0 = A*tAm[na] + tCm[na], ctp1 = A*tAp[nb] + tCp[nb];
        const int o = l*ndir;
        const bool lo = (na>0    && -neg(ctm0)*rdth*double(Nl[na])>thr[l]*tdi[o+na]);
        const bool hi = (nb<nq-1 &&  pos(ctp1)*rdth*double(Nl[nb])>thr[l]*tdi[o+nb]);
        if(!lo && !hi)
        return;

        const int w = std::max(2,nb-na+1);
        const int na2 = lo ? std::max(0,na-w) : na, nb2 = hi ? std::min(nq-1,nb+w) : nb;

            if(src!=nullptr)
            {
            src->set_rows(rsum.data());
            src->compute(N,d,kc,cgc,P.data(),D.data(),ma+na2,ma+nb2,l,l);
            }

        wlo[l] = na2;
        whi[l] = nb2;
        wmin = std::min(wmin,na2);
        wmax = std::max(wmax,nb2);
        assemble(l,na2,nb2);
        solve1(l);
        store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
        }
    };

    // C: forward and back substitution. Without frequency shift the frequencies are independent;
    // with it, the inflow from l-1 takes the value just solved (Gauss-Seidel in sigma)
    if(sparse)
    for(int l=lmin; l<=lmax; ++l)
    {
    solve1(l);
    store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
    widen(l);
    }
    else if(!fs)
    {
        for(int n=wmin; n<=wmax; ++n)
        for(int l=lmin; l<=lmax; ++l)
        {
        if(n<wlo[l] || n>whi[l])
        continue;

        const int x = l*ndir + n;
        tdp[x] = (trh[x] - (n>wlo[l] ? tla[x]*tdp[x-1] : 0.0))*trd[x];
        }

        for(int n=wmax; n>=wmin; --n)
        for(int l=lmin; l<=lmax; ++l)
        {
        if(n<wlo[l] || n>whi[l])
        continue;

        const int x = l*ndir + n;
        tso[x] = (n<whi[l]) ? tdp[x] - tcp[x]*tso[x+1] : tdp[x];
        }

        for(int l=lmin; l<=lmax; ++l)
        store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
    }
    else
    for(int l=lmin; l<=lmax; ++l)
    {
    const int o = l*ndir, na = wlo[l], nb = whi[l];

        if(na<=nb)
        {
            if(l>0)
            {
            const float *Nlm = N + g.bin(l-1,ma);
            for(int n=na; n<=nb; ++n)
            trh[o+n] += tsi[o+n]*double(Nlm[n]);
            }

            for(int n=na; n<=nb; ++n)
            tdp[o+n] = (trh[o+n] - (n>na ? tla[o+n]*tdp[o+n-1] : 0.0))*trd[o+n];

            tso[o+nb] = tdp[o+nb];
            for(int n=nb-1; n>=na; --n)
            tso[o+n] = tdp[o+n] - tcp[o+n]*tso[o+n+1];
        }

    store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
    }

    // new active ranges of the cell
    if(sparse)
    keep(p,e,q,ic,jc);

    i = ic;
    j = jc;
}

/*--------------------------------------------------------------------
Surfbeat (A 775 2): one frequency, Crank-Nicolson in time for the
transport (theta = 1/2: no numerical diffusion from the time stepping,
the wave groups keep their shape over long distances), second-order
geographic fluxes (van Leer, deferred correction), the source terms
implicit (Patankar). Old values: N0 (the spectra at the start of the
step, halo included) and the boundary rows of the previous step.
Positive for c dt/dx <= 2 (the host's time step is limited by the
long-wave celerity, which exceeds c_g); N is clipped at zero.
--------------------------------------------------------------------*/

void seastate_implicit::cell_surfbeat(lexer *p, fdm_seastate *e, int q, const seastate_store *N0s, double rdt,
                                      const vector<float> &Nb, const int side[4], bool refraction)
{
    const seastate_grid &g = *e->grid;
    const double th = 0.5, th1 = 1.0-th;
    const int nbin = g.nbin;

    float *N = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cgc = e->cg->spec(i,j);

    const int ic = i, jc = j;
    const float *Nc0 = (N0s!=nullptr) ? N0s->spec(ic,jc) : N;

    neighbour W,E,S,Nn;
    set_neighbour(p,e,ic-1,jc,Nb,Nbs,side,W);
    set_neighbour(p,e,ic+1,jc,Nb,Nbs,side,E);
    set_neighbour(p,e,ic,jc-1,Nb,Nbs,side,S);
    set_neighbour(p,e,ic,jc+1,Nb,Nbs,side,Nn);

    // old values of the neighbours: active cells from N0, the x- rows of the previous step, Nb
    auto old = [&](const neighbour &nb, int ni, int nj) -> const float*
    {
        if(nb.N==nullptr)
        return nullptr;

        if(nb.cg!=nullptr)
        return (N0s!=nullptr) ? N0s->spec(ni,nj) : nb.N;

        if(ni+p->origin_i<0 && Nbs0[0]!=nullptr && nj>=0 && nj<p->knoy)
        return Nbs0[0]->data() + size_t(nj)*size_t(nbin);

        return nb.N;
    };

    const float *W0 = old(W,ic-1,jc), *E0 = old(E,ic+1,jc), *S0 = old(S,ic,jc-1), *N0n = old(Nn,ic,jc+1);

    // second order: active neighbours and active cells two cells away, new and old values
    const float *aW = W.cg ? W.N : nullptr, *aE = E.cg ? E.N : nullptr, *aS = S.cg ? S.N : nullptr, *aN = Nn.cg ? Nn.N : nullptr;
    const float *oW = W.cg ? W0 : nullptr, *oE = E.cg ? E0 : nullptr, *oS = S.cg ? S0 : nullptr, *oN = Nn.cg ? N0n : nullptr;
    const float *aWW = active(e,ic-2,jc), *aEE = active(e,ic+2,jc), *aSS = active(e,ic,jc-2), *aNN = active(e,ic,jc+2);
    const float *oWW = (aWW && N0s) ? N0s->spec(ic-2,jc) : aWW, *oEE = (aEE && N0s) ? N0s->spec(ic+2,jc) : aEE;
    const float *oSS = (aSS && N0s) ? N0s->spec(ic,jc-2) : aSS, *oNN = (aNN && N0s) ? N0s->spec(ic,jc+2) : aNN;

    const double rdx = 1.0/p->DXN[IP];
    const double rdy = 1.0/p->DYN[JP];
    const double rdth = 1.0/g.dtheta;

    const double d = e->depth(ic,jc);
    const double U = e->U(ic,jc), V = e->V(ic,jc);
    const bool kin = (e->refr(ic,jc)==1);

    const double ddx = e->ddx(ic,jc), ddy = e->ddy(ic,jc);
    const double dUdx = e->dUdx(ic,jc), dUdy = e->dUdy(ic,jc), dVdx = e->dVdx(ic,jc), dVdy = e->dVdy(ic,jc);

    const int ma = m0[q], mb = m1[q], nq = mb-ma+1;

    if(src!=nullptr)
    src->compute(N,d,kc,cgc,P.data(),D.data());

    const double sig = g.sig[0];
    const double cgl = cgc[0];
    const double A = kin ? seastate_refraction(sig,kc[0],d) : 0.0;

    for(int n=0; n<nq; ++n)
    {
    const int m = ma+n;
    const int b = g.bin(0,m);
    const double cs = g.costh[m], sn = g.sinth[m];
    const int mm = (m-1+ndir)%ndir, mp = (m+1)%ndir;

    // c_theta at the faces m+1/2 and m-1/2
    double ctp=0.0, ctm=0.0;

        if(refraction && kin)
        {
        const double sp = g.sinthf[m],  cp_ = g.costhf[m];
        const double sm = g.sinthf[mm], cm  = g.costhf[mm];

        ctp = A*(sp*ddx - cp_*ddy) + sp*cp_*(dUdx-dVdy) + sp*sp*dVdx - cp_*cp_*dUdy;
        ctm = A*(sm*ddx - cm*ddy)  + sm*cm*(dUdx-dVdy)  + sm*sm*dVdx - cm*cm*dUdy;
        }

    // geographic face velocities
    const double cxc = cgl*cs + U;
    const double cyc = cgl*sn + V;

    const double cxw = W.cg  ? 0.5*(cxc + double(W.cg[0])*cs + W.U)   : cxc;
    const double cxe = E.cg  ? 0.5*(cxc + double(E.cg[0])*cs + E.U)   : cxc;
    const double cys = S.cg  ? 0.5*(cyc + double(S.cg[0])*sn + S.V)   : cyc;
    const double cyn = Nn.cg ? 0.5*(cyc + double(Nn.cg[0])*sn + Nn.V) : cyc;

    const double out  = (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy;
    const double outt = (pos(ctp) - neg(ctm))*rdth;

    di[n] = rdt + th*(out + outt) + (src!=nullptr ? D[b] : 0.0);

    double r = rdt*double(Nc0[b]) - th1*(out + outt)*double(Nc0[b]);

    if(src!=nullptr)
    r += P[b];

    // zero-gradient sides: inflow of the cell's own spectrum
    double aself = 0.0;
    if(W.self)  aself += pos(cxw)*rdx;
    if(E.self)  aself -= neg(cxe)*rdx;
    if(S.self)  aself += pos(cys)*rdy;
    if(Nn.self) aself -= neg(cyn)*rdy;

        if(aself>0.0)
        {
        r += th1*aself*double(Nc0[b]);

            if(di[n]-th*aself>=0.5*di[n])
            di[n] -= th*aself;
            else
            r += th*aself*double(N[b]);
        }

    // inflow from the neighbours, new (latest) and old values
    if(W.N)  r += pos(cxw)*rdx*(th*W.N[b]  + th1*W0[b]);
    if(E.N)  r -= neg(cxe)*rdx*(th*E.N[b]  + th1*E0[b]);
    if(S.N)  r += pos(cys)*rdy*(th*S.N[b]  + th1*S0[b]);
    if(Nn.N) r -= neg(cyn)*rdy*(th*Nn.N[b] + th1*N0n[b]);

    // second-order correction of the geographic fluxes
    {
    const double fe = cxe>0.0 ? tvd(cxe,aW,N,aE,b)  : tvd(cxe,aEE,aE,N,b);
    const double fw = cxw>0.0 ? tvd(cxw,aWW,aW,N,b) : tvd(cxw,aE,N,aW,b);
    const double fn = cyn>0.0 ? tvd(cyn,aS,N,aN,b)  : tvd(cyn,aNN,aN,N,b);
    const double fs = cys>0.0 ? tvd(cys,aSS,aS,N,b) : tvd(cys,aN,N,aS,b);

    const double ge = cxe>0.0 ? tvd(cxe,oW,Nc0,oE,b)  : tvd(cxe,oEE,oE,Nc0,b);
    const double gw = cxw>0.0 ? tvd(cxw,oWW,oW,Nc0,b) : tvd(cxw,oE,Nc0,oW,b);
    const double gn = cyn>0.0 ? tvd(cyn,oS,Nc0,oN,b)  : tvd(cyn,oNN,oN,Nc0,b);
    const double gs = cys>0.0 ? tvd(cys,oSS,oS,Nc0,b) : tvd(cys,oN,Nc0,oS,b);

    r -= th*((fe - fw)*rdx + (fn - fs)*rdy) + th1*((ge - gw)*rdx + (gn - gs)*rdy);
    }

    // theta: within the quadrant implicit (tridiagonal), outside with the latest values; old values explicit
    r += th1*(pos(ctm)*double(Nc0[g.bin(0,mm)]) - neg(ctp)*double(Nc0[g.bin(0,mp)]))*rdth;

    la[n] = up[n] = 0.0;

        if(n>0)
        la[n] = -th*pos(ctm)*rdth;
        else
        r += th*pos(ctm)*rdth*N[g.bin(0,mm)];

        if(n<nq-1)
        up[n] = th*neg(ctp)*rdth;
        else
        r -= th*neg(ctp)*rdth*N[g.bin(0,mp)];

    rhs[n] = r;
    }

    // Thomas algorithm
    for(int n=0; n<nq; ++n)
    {
    const double den = di[n] - (n>0 ? la[n]*cp[n-1] : 0.0);

        if(!(den>1.0e-300))
        {
        for(int nn=0; nn<nq; ++nn)
        N[g.bin(0,ma+nn)] = 0.0f;

        i = ic;
        j = jc;
        return;
        }

    cp[n] = up[n]/den;
    dp[n] = (rhs[n] - (n>0 ? la[n]*dp[n-1] : 0.0))/den;
    }

    sol[nq-1] = dp[nq-1];
    for(int n=nq-2; n>=0; --n)
    sol[n] = dp[n] - cp[n]*sol[n+1];

    for(int n=0; n<nq; ++n)
    N[g.bin(0,ma+n)] = float(std::max(sol[n],0.0));

    i = ic;
    j = jc;
}
