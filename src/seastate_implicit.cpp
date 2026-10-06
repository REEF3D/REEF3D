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
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

seastate_implicit::seastate_implicit(lexer *p, fdm_seastate *e) : src(nullptr), Nbx(nullptr), Nbx0(nullptr), second(false)
{
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

    csig.assign(nsig,0.0);
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

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
    i = idown ? p->knox-1-ii : ii;
    j = jdown ? p->knoy-1-jj : jj;

        if(e->wet(i,j)==1)
        {
            if(second)
            cell_surfbeat(p,e,q,N0,rdt,Nb,side,refraction);
            else
            cell(p,e,q,N0!=nullptr ? N0->spec(i,j) : nullptr,rdt,Nb,side,refraction,fshift);
        }
    }
}

namespace
{
    // neighbour of a cell: active cell, boundary with incoming spectrum, zero-gradient side, or nothing (land, open side)
    struct neighbour
    {
        const float *N = nullptr;     // spectrum (cell or boundary), nullptr: no inflow
        const float *cg = nullptr;    // group velocity, nullptr: no active cell
        double U = 0.0, V = 0.0;
        bool self = false;            // zero-gradient side: inflow of the cell's own spectrum
    };

    void set_neighbour(lexer *p, fdm_seastate *e, int ni, int nj, const vector<float> &Nb, const vector<float> *Nbx, const int side[4], neighbour &nb)
    {
        nb = neighbour();

        if(e->wet(ni,nj)==1)
        {
        nb.N  = e->N->spec(ni,nj);
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

        // surfbeat: x- boundary spectrum of the row
        if(s==0 && side[0]==1 && Nbx!=nullptr && nj>=0 && nj<p->knoy)
        nb.N = Nbx->data() + size_t(nj)*size_t(e->N->nbins());

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

void seastate_implicit::cell(lexer *p, fdm_seastate *e, int q, const float *N0, double rdt,
                             const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    const seastate_grid &g = *e->grid;

    float *N = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cgc = e->cg->spec(i,j);

    const int ic = i, jc = j;

    neighbour W,E,S,Nn;
    set_neighbour(p,e,ic-1,jc,Nb,Nbx,side,W);
    set_neighbour(p,e,ic+1,jc,Nb,Nbx,side,E);
    set_neighbour(p,e,ic,jc-1,Nb,Nbx,side,S);
    set_neighbour(p,e,ic,jc+1,Nb,Nbx,side,Nn);

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


    // source terms from the latest spectrum of the cell
    if(src!=nullptr)
    src->compute(N,d,kc,cgc,P.data(),D.data());

    const double *lim = (src!=nullptr) ? src->limit() : nullptr;

    for(int l=0; l<nsig; ++l)
    {
    const double sig = g.sig[l];
    const double cgl = cgc[l];
    const double A = kin ? seastate_refraction(sig,kc[l],d) : 0.0;
    const double rdsig = 1.0/g.dsig[l];

        for(int n=0; n<nq; ++n)
        {
        const int m = ma+n;
        const int b = g.bin(l,m);
        const double cs = g.costh[m], sn = g.sinth[m];

        // c_sigma at l-1, l, l+1 for this direction
        double csm=0.0, csc=0.0, csp=0.0;

            if(fshift && kin)
            {
            const double cur = cs*cs*dUdx + sn*cs*(dUdy+dVdx) + sn*sn*dVdy;
            const double dep = dddt + U*ddx + V*ddy;

            for(int ll=std::max(l-1,0); ll<=std::min(l+1,nsig-1); ++ll)
            {
            const double Al = seastate_refraction(g.sig[ll],kc[ll],d);
            const double c = kc[ll]*Al*dep - double(cgc[ll])*kc[ll]*cur;

            if(ll==l-1) csm=c;
            if(ll==l)   csc=c;
            if(ll==l+1) csp=c;
            }
            }

        const double csu = (l<nsig-1) ? 0.5*(csc+csp) : csc;     // face l+1/2
        const double csl = (l>0)      ? 0.5*(csc+csm) : csc;     // face l-1/2

        // c_theta at the faces m+1/2 and m-1/2
        double ctp=0.0, ctm=0.0;

            if(refraction && kin)
            {
            const int mm = (m-1+ndir)%ndir;

            const double sp = g.sinthf[m],  cp_ = g.costhf[m];
            const double sm = g.sinthf[mm], cm  = g.costhf[mm];

            ctp = A*(sp*ddx - cp_*ddy) + sp*cp_*(dUdx-dVdy) + sp*sp*dVdx - cp_*cp_*dUdy;
            ctm = A*(sm*ddx - cm*ddy)  + sm*cm*(dUdx-dVdy)  + sm*sm*dVdx - cm*cm*dUdy;
            }

        // geographic face velocities
        const double cxc = cgl*cs + U;
        const double cyc = cgl*sn + V;

        const double cxw = W.cg  ? 0.5*(cxc + double(W.cg[l])*cs + W.U)   : cxc;
        const double cxe = E.cg  ? 0.5*(cxc + double(E.cg[l])*cs + E.U)   : cxc;
        const double cys = S.cg  ? 0.5*(cyc + double(S.cg[l])*sn + S.V)   : cyc;
        const double cyn = Nn.cg ? 0.5*(cyc + double(Nn.cg[l])*sn + Nn.V) : cyc;

        // diagonal: 1/dt + outflow (+ D)
        const double dg = rdt + (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy
                              + (pos(csu) - neg(csl))*rdsig + (src!=nullptr ? D[b] : 0.0);

        di[n] = dg + (pos(ctp) - neg(ctm))*rdth;

        // right-hand side: old time level + inflow from neighbours (+ P)
        double r = rdt*(N0!=nullptr ? double(N0[b]) : double(N[b]));

        if(src!=nullptr)
        r += P[b];

        // zero-gradient sides: inflow of the cell's own spectrum, implicit while the
        // non-theta part of the diagonal stays at least half of its value (M-matrix),
        // otherwise with the latest value
        double aself = 0.0;
        if(W.self)  aself += pos(cxw)*rdx;
        if(E.self)  aself -= neg(cxe)*rdx;
        if(S.self)  aself += pos(cys)*rdy;
        if(Nn.self) aself -= neg(cyn)*rdy;

            if(aself>0.0)
            {
                if(dg-aself>=0.5*dg)
                di[n] -= aself;
                else
                r += aself*double(N[b]);
            }

        if(W.N)  r += pos(cxw)*rdx*W.N[b];
        if(E.N)  r -= neg(cxe)*rdx*E.N[b];
        if(S.N)  r += pos(cys)*rdy*S.N[b];
        if(Nn.N) r -= neg(cyn)*rdy*Nn.N[b];

        if(l>0)
        r += pos(csl)*rdsig*N[g.bin(l-1,m)];

        if(l<nsig-1)
        r -= neg(csu)*rdsig*N[g.bin(l+1,m)];

        // theta: within the quadrant implicit (tridiagonal), outside with the latest values
        la[n] = up[n] = 0.0;

            if(n>0)
            la[n] = -pos(ctm)*rdth;
            else
            r += pos(ctm)*rdth*N[g.bin(l,(m-1+ndir)%ndir)];

            if(n<nq-1)
            up[n] = neg(ctp)*rdth;
            else
            r -= neg(ctp)*rdth*N[g.bin(l,(m+1)%ndir)];

        rhs[n] = r;
        }

        // Thomas algorithm (the tridiagonal block of an M-matrix: no pivoting needed)
        bool ok = true;

        for(int n=0; n<nq; ++n)
        {
        const double den = di[n] - (n>0 ? la[n]*cp[n-1] : 0.0);

            if(!(den>1.0e-300))
            {
            ok = false;
            break;
            }

        cp[n] = up[n]/den;
        dp[n] = (rhs[n] - (n>0 ? la[n]*dp[n-1] : 0.0))/den;
        }

        if(!ok)
        {
        for(int n=0; n<nq; ++n)
        N[g.bin(l,ma+n)] = 0.0f;

        continue;
        }

        sol[nq-1] = dp[nq-1];
        for(int n=nq-2; n>=0; --n)
        sol[n] = dp[n] - cp[n]*sol[n+1];

        // action density limiter (wind sea, seastate_source): |change per iteration| <= dNmax
        if(lim!=nullptr)
        for(int n=0; n<nq; ++n)
        {
        const double o = double(N[g.bin(l,ma+n)]);
        sol[n] = std::min(std::max(sol[n],o-lim[l]),o+lim[l]);
        }

        for(int n=0; n<nq; ++n)
        N[g.bin(l,ma+n)] = float(sol[n]);
    }

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
    set_neighbour(p,e,ic-1,jc,Nb,Nbx,side,W);
    set_neighbour(p,e,ic+1,jc,Nb,Nbx,side,E);
    set_neighbour(p,e,ic,jc-1,Nb,Nbx,side,S);
    set_neighbour(p,e,ic,jc+1,Nb,Nbx,side,Nn);

    // old values of the neighbours: active cells from N0, the x- rows of the previous step, Nb
    auto old = [&](const neighbour &nb, int ni, int nj) -> const float*
    {
        if(nb.N==nullptr)
        return nullptr;

        if(nb.cg!=nullptr)
        return (N0s!=nullptr) ? N0s->spec(ni,nj) : nb.N;

        if(ni+p->origin_i<0 && Nbx0!=nullptr && nj>=0 && nj<p->knoy)
        return Nbx0->data() + size_t(nj)*size_t(nbin);

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
