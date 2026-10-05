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

#include"spectral_implicit.h"
#include"spectral_exchange.h"
#include"spectral_grid.h"
#include"spectral_store.h"
#include"spectral_dispersion.h"
#include"fdm_spectral.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

spectral_implicit::spectral_implicit(lexer *p, fdm_spectral *e)
{
    const spectral_grid &g = *e->grid;

    nsig = g.nsig;
    ndir = g.ndir;

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

void spectral_implicit::iterate(lexer *p, ghostcell *pgc, fdm_spectral *e, spectral_exchange *pex, const spectral_store *N0,
                                double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    for(int q=0; q<4; ++q)
    {
    sweep(p,e,q,N0,rdt,Nb,side,refraction,fshift);

    if(pex!=nullptr)
    pex->start(p,pgc,*e->N);
    }
}

void spectral_implicit::sweep(lexer *p, fdm_spectral *e, int q, const spectral_store *N0, double rdt,
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
        cell(p,e,q,N0!=nullptr ? N0->spec(i,j) : nullptr,rdt,Nb,side,refraction,fshift);
    }
}

namespace
{
    // neighbour of a cell: active cell, boundary with incoming spectrum, or nothing (land, open side)
    struct neighbour
    {
        const float *N = nullptr;     // spectrum (cell or boundary), nullptr: no inflow
        const float *cg = nullptr;    // group velocity, nullptr: no active cell
        double U = 0.0, V = 0.0;
    };

    void set_neighbour(lexer *p, fdm_spectral *e, int ni, int nj, const vector<float> &Nb, const int side[4], neighbour &nb)
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
    }

    inline double pos(double a) {return a>0.0 ? a : 0.0;}
    inline double neg(double a) {return a<0.0 ? a : 0.0;}
}

void spectral_implicit::cell(lexer *p, fdm_spectral *e, int q, const float *N0, double rdt,
                             const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    const spectral_grid &g = *e->grid;

    float *N = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cgc = e->cg->spec(i,j);

    const int ic = i, jc = j;

    neighbour W,E,S,Nn;
    set_neighbour(p,e,ic-1,jc,Nb,side,W);
    set_neighbour(p,e,ic+1,jc,Nb,side,E);
    set_neighbour(p,e,ic,jc-1,Nb,side,S);
    set_neighbour(p,e,ic,jc+1,Nb,side,Nn);

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

    for(int l=0; l<nsig; ++l)
    {
    const double sig = g.sig[l];
    const double cgl = cgc[l];
    const double A = kin ? spectral_refraction(sig,kc[l],d) : 0.0;
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
            const double Al = spectral_refraction(g.sig[ll],kc[ll],d);
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

        // diagonal: 1/dt + outflow
        di[n] = rdt + (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy
                    + (pos(csu) - neg(csl))*rdsig + (pos(ctp) - neg(ctm))*rdth;

        // right-hand side: old time level + inflow from neighbours
        double r = rdt*(N0!=nullptr ? double(N0[b]) : double(N[b]));

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

        for(int n=0; n<nq; ++n)
        N[g.bin(l,ma+n)] = float(sol[n]);
    }

    i = ic;
    j = jc;
}
