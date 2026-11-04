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

#include"seastate_f.h"
#include"seastate_exchange.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_implicit.h"
#include"slice4.h"
#include"seastate_store.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

/*--------------------------------------------------------------------
GSE closure (A 776, Phase 25), split from the transport:
after the iterations of a nonstationary step, the flux inside each bin
as a diffusion over the step,
  (1 - dt div D grad) N = N*,   D = t_a Sigma,
  Sigma = diag(dcg^2, (cg dtheta)^2)/12 along and across the direction
  of the bin (the spread of its velocity), t_a the age (A 776 1: A 777,
  A 776 2: the model time),
implicit in the axis terms (A 778 red-black Gauss-Seidel iterations
over the cells by their global index, the halo exchanged before each
colour: the same result on any decomposition; the matrix diagonally
dominant, N stays positive), the cross term explicit from N* (the new N
clipped at zero), the faces to dry cells closed. Validation 79: a swell
field after 3 days against a continuous spectrum, L2 0.46-0.57 without
the closure (second order), 0.136 with 2 iterations (+34 % time), 0.129
with 8. A 776 3 (Phase 26): the age per bin, transported (age_step), the
faces with the energy-weighted age of the two cells. The coupled form in the sweeps
needs the diffusion converged within the step (validation 79: one
iteration per step lags, the old time level slows the propagation).
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
The transported age (A 776 3, Phase 26): M = t_a N per bin, with the
transport of the spectra (the solver of the step without source terms,
the same refraction, frequency shift, obstacles and pipeline):
  dM/dt + div(c M) + ... = N   (every bin ages with the model time)
from M + dt N at the start of the step; the boundary and the initial
spectra with the age A 777. The source terms are not in it: a gain of
N makes the bin younger (M/N falls), a loss older; so in the wind sea,
where 1.7 U10 cos(theta - theta_w) > c (the wind-sea criterion of
WAVEWATCH III), M = 0 (the closure is for the swell), and t_a is at most
the model time plus A 777. The GSE step diffuses M with the same
coefficients as N.
--------------------------------------------------------------------*/

void seastate_f::age_ini(lexer *p, ghostcell *pgc)
{
    const int nbin = e->grid->nbin;
    auto mk = [&]() {seastate_store *s = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,nbin,p->A704); s->build(e->wet.V); return s;};
    gM = mk();
    gM0 = mk();
    gMs = mk();

    pasolv = new seastate_implicit(p,e);
    pasolv->geographic_order(p->A796);
    pasolv->threads(p->A798);
    if(bser!=nullptr)
    {
    const std::vector<float> *s[4] = {&Mside[0],&Mside[1],&Mside[2],&Mside[3]};
    pasolv->boundary_sides(s);
    }
    if(dca!=nullptr)
    pasolv->diffraction(dca,dcax,dcay);
    if(pobs!=nullptr)
    pasolv->obstacles(pobs);
    if(kba)
    pasolv->pipeline(p,pgc);
}

// every spectrum of a store: a = s*x + y (y may be nullptr)
static void axpy(lexer *p, seastate_store &a, double s, const seastate_store &x, const seastate_store *y, int nbin)
{
    for(int i=p->imin; i<p->imin+p->imax; ++i)
    for(int j=p->jmin; j<p->jmin+p->jmax; ++j)
    {
    float *A = a.spec(i,j);
    const float *X = x.spec(i,j);
    if(A==nullptr || X==nullptr)
    continue;
    const float *Y = (y!=nullptr) ? y->spec(i,j) : nullptr;
        for(int b=0; b<nbin; ++b)
        A[b] = float(s*double(X[b]) + (Y!=nullptr ? double(Y[b]) : 0.0));
    }
}

void seastate_f::age_start(lexer *p)
{
    const int nbin = e->grid->nbin;
    if(!gMset)
    {
    axpy(p,*gM,p->A777,*e->N,nullptr,nbin);
    gMset = true;
    }
    // the old time level with the ageing: M + dt N
    axpy(p,*gM0,dtw,*e->N,gM,nbin);
}

void seastate_f::age_step(lexer *p, ghostcell *pgc)
{
    const seastate_grid &g = *e->grid;
    const int nsig = g.nsig, ndir = g.ndir, nbin = g.nbin;
    const double rdt = 1.0/dtw;
    const bool refraction = (p->A713==1);
    const bool fshift = (p->A714==1 && sb==nullptr);

    // the boundary spectra with their age
    Mb.resize(Nb.size());
    for(size_t k=0; k<Nb.size(); ++k)
    Mb[k] = float(p->A777*double(Nb[k]));
    for(int k=0; k<4; ++k)
    {
    Mside[k].resize(Nside[k].size());
    for(size_t n=0; n<Nside[k].size(); ++n)
    Mside[k][n] = float(p->A777*double(Nside[k][n]));
    }

    // the transport of M: the solver of the step on M (the spectra swapped), one iteration
    seastate_store *Nk = e->N;
    e->N = gM;
    pasolv->iterate(p,pgc,e,pex,gM0,rdt,Mb,side,refraction,fshift);
    e->N = Nk;

    // the wind sea: no age (the closure is for the swell); the age at most the model time plus A 777
    const double tmax = p->simtime + dtw + p->A777;
    const bool wind = (p->A730>0);
        for(int i=p->imin; i<p->imin+p->imax; ++i)
        for(int j=p->jmin; j<p->jmin+p->jmax; ++j)
        {
        float *M = gM->spec(i,j);
        const float *N = e->N->spec(i,j);
        if(M==nullptr || N==nullptr)
        continue;
        const bool own = (i>=0 && j>=0 && i<p->knox && j<p->knoy);
        double U = 0.0, cw = 1.0, sw = 0.0;
            if(wind && own)
            {
            const double wd = (wU10!=nullptr) ? double((*wdir)(i,j)) : p->A731_dir*3.14159265358979323846/180.0;
            U = (wU10!=nullptr) ? double((*wU10)(i,j)) : p->A731_u10;
            cw = std::cos(wd); sw = std::sin(wd);
            }
        const float *kc = e->kw->spec(i,j);
            for(int l=0; l<nsig; ++l)
            {
            const double c = (kc!=nullptr && kc[l]>0.0f) ? g.sig[l]/double(kc[l]) : 1.0e30;
                for(int m=0; m<ndir; ++m)
                {
                const int b = g.bin(l,m);
                double v = std::max(double(M[b]),0.0);
                if(U>0.0 && 1.7*U*(g.costh[m]*cw + g.sinth[m]*sw)>c)
                v = 0.0;
                v = std::min(v,tmax*std::max(double(N[b]),0.0));
                M[b] = float(v);
                }
            }
        }
    pex->start(p,pgc,*gM);
}

void seastate_f::gse_step(lexer *p, ghostcell *pgc)
{
    if(p->A776==0 || p->A700!=1)
    return;

    const double dt = dtw;
    if(p->A776==3)
    {
    gse_age_step(p,pgc);
    return;
    }
    const double ta = (p->A776==2) ? p->simtime + dt : p->A777;
    if(!(ta>0.0) || !(dt>0.0))
    return;

    const seastate_grid &g = *e->grid;
    const int nsig = g.nsig, ndir = g.ndir, nbin = g.nbin;

    if(gNs==nullptr)
    {
    gNs = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,nbin,p->A704);
    gNs->build(e->wet.V);
    }
    gNs->copy_from(*e->N);

    const int nit = p->A778;
    auto wet = [&](int ii, int jj) {return e->wet(ii,jj)==1 && e->N->spec(ii,jj)!=nullptr;};

    // red-black Gauss-Seidel by the global index of the cell: the same result on any decomposition
    for(int it=0; it<nit; ++it)
    for(int colour=0; colour<2; ++colour)
    {
    pex->start(p,pgc,*e->N);

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        {
            if(!wet(i,j) || ((i+p->origin_i+j+p->origin_j)&1)!=colour)
            continue;

        float *N = e->N->spec(i,j);
        const float *Ns = gNs->spec(i,j);
        const float *cg = e->cg->spec(i,j);
        const double rdx = 1.0/p->DXN[i+marge], rdy = (p->j_dir==1) ? 1.0/p->DYN[j+marge] : 0.0;
        const double rx2 = rdx*rdx, ry2 = rdy*rdy, rxy = 0.25*rdx*rdy;

        const float *nW = wet(i-1,j) ? e->N->spec(i-1,j) : nullptr;
        const float *nE = wet(i+1,j) ? e->N->spec(i+1,j) : nullptr;
        const float *nS = (p->j_dir==1 && wet(i,j-1)) ? e->N->spec(i,j-1) : nullptr;
        const float *nN = (p->j_dir==1 && wet(i,j+1)) ? e->N->spec(i,j+1) : nullptr;
        // the cross term from N* of the corners (all four wet)
        const bool cr = (p->j_dir==1 && wet(i+1,j+1) && wet(i-1,j+1) && wet(i+1,j-1) && wet(i-1,j-1));
        const float *sNE = cr ? gNs->spec(i+1,j+1) : nullptr, *sNW = cr ? gNs->spec(i-1,j+1) : nullptr;
        const float *sSE = cr ? gNs->spec(i+1,j-1) : nullptr, *sSW = cr ? gNs->spec(i-1,j-1) : nullptr;

            for(int l=0; l<nsig; ++l)
            {
            const int la = std::max(l-1,0), lb = std::min(l+1,nsig-1);
            const double dcg = (lb>la) ? std::fabs(double(cg[lb]) - double(cg[la]))/(g.sig[lb] - g.sig[la])*g.dsig[l] : 0.0;
            const double css = ta*dcg*dcg/12.0, cnn = ta*std::pow(double(cg[l])*g.dtheta,2.0)/12.0;

                for(int m=0; m<ndir; ++m)
                {
                const int b = g.bin(l,m);
                const double cs = g.costh[m], sn = g.sinth[m];
                const double dxx = css*cs*cs + cnn*sn*sn, dyy = css*sn*sn + cnn*cs*cs, dxy = (css-cnn)*sn*cs;
                double den = 1.0, num = double(Ns[b]);
                    if(nW) {den += dt*dxx*rx2; num += dt*dxx*rx2*double(nW[b]);}
                    if(nE) {den += dt*dxx*rx2; num += dt*dxx*rx2*double(nE[b]);}
                    if(nS) {den += dt*dyy*ry2; num += dt*dyy*ry2*double(nS[b]);}
                    if(nN) {den += dt*dyy*ry2; num += dt*dyy*ry2*double(nN[b]);}
                    if(cr)
                    num += dt*2.0*dxy*rxy*(double(sNE[b]) - double(sNW[b]) - double(sSE[b]) + double(sSW[b]));
                N[b] = float(std::max(num,0.0)/den);
                }
            }
        }
    }

    pex->start(p,pgc,*e->N);
}

// A 776 3 (Phase 26): the closure with the transported age; every face with the energy-weighted age of its two cells
// (M_i + M_n)/(N_i + N_n) and the mean Sigma of the two (the face coefficients symmetric: the step conserves action);
// M diffused with the same coefficients; the cross term with the age of the cell
void seastate_f::gse_age_step(lexer *p, ghostcell *pgc)
{
    const double dt = dtw;
    if(!(dt>0.0))
    return;
    const seastate_grid &g = *e->grid;
    const int nsig = g.nsig, ndir = g.ndir, nbin = g.nbin;

    if(gNs==nullptr)
    {
    gNs = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,nbin,p->A704);
    gNs->build(e->wet.V);
    }
    gNs->copy_from(*e->N);
    gMs->copy_from(*gM);

    auto wet = [&](int ii, int jj) {return e->wet(ii,jj)==1 && e->N->spec(ii,jj)!=nullptr;};
    // Sigma/t (along, across) of a cell and a frequency
    auto sig2 = [&](const float *cg, int l, double &ss, double &nn)
    {
    const int la = std::max(l-1,0), lb = std::min(l+1,nsig-1);
    const double dcg = (lb>la) ? std::fabs(double(cg[lb]) - double(cg[la]))/(g.sig[lb] - g.sig[la])*g.dsig[l] : 0.0;
    ss = dcg*dcg/12.0; nn = std::pow(double(cg[l])*g.dtheta,2.0)/12.0;
    };
    std::vector<double> sW(nsig*2), sE(nsig*2), sS(nsig*2), sN(nsig*2), sC(nsig*2);

    for(int it=0; it<p->A778; ++it)
    for(int colour=0; colour<2; ++colour)
    {
    pex->start(p,pgc,*e->N);
    pex->start(p,pgc,*gM);

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        {
            if(!wet(i,j) || ((i+p->origin_i+j+p->origin_j)&1)!=colour)
            continue;

        float *N = e->N->spec(i,j), *M = gM->spec(i,j);
        const float *Ns = gNs->spec(i,j), *Ms = gMs->spec(i,j);
        const double rdx = 1.0/p->DXN[i+marge], rdy = (p->j_dir==1) ? 1.0/p->DYN[j+marge] : 0.0;
        const double rx2 = rdx*rdx, ry2 = rdy*rdy, rxy = 0.25*rdx*rdy;

        struct nb {const float *N, *M, *Ns, *Ms, *cg;};
        auto take = [&](int ii, int jj, bool ok) {nb r = {nullptr,nullptr,nullptr,nullptr,nullptr};
            if(ok && wet(ii,jj)) {r.N = e->N->spec(ii,jj); r.M = gM->spec(ii,jj); r.Ns = gNs->spec(ii,jj); r.Ms = gMs->spec(ii,jj); r.cg = e->cg->spec(ii,jj);}
            return r;};
        const nb W = take(i-1,j,true), E = take(i+1,j,true), S = take(i,j-1,p->j_dir==1), Nn = take(i,j+1,p->j_dir==1);
        const bool cr = (p->j_dir==1 && wet(i+1,j+1) && wet(i-1,j+1) && wet(i+1,j-1) && wet(i-1,j-1));
        const float *sNE = cr ? gNs->spec(i+1,j+1) : nullptr, *sNW = cr ? gNs->spec(i-1,j+1) : nullptr;
        const float *sSE = cr ? gNs->spec(i+1,j-1) : nullptr, *sSW = cr ? gNs->spec(i-1,j-1) : nullptr;
        const float *mNE = cr ? gMs->spec(i+1,j+1) : nullptr, *mNW = cr ? gMs->spec(i-1,j+1) : nullptr;
        const float *mSE = cr ? gMs->spec(i+1,j-1) : nullptr, *mSW = cr ? gMs->spec(i-1,j-1) : nullptr;

        const float *cgc = e->cg->spec(i,j);
            for(int l=0; l<nsig; ++l)
            {
            sig2(cgc,l,sC[2*l],sC[2*l+1]);
            if(W.cg)  sig2(W.cg,l,sW[2*l],sW[2*l+1]);
            if(E.cg)  sig2(E.cg,l,sE[2*l],sE[2*l+1]);
            if(S.cg)  sig2(S.cg,l,sS[2*l],sS[2*l+1]);
            if(Nn.cg) sig2(Nn.cg,l,sN[2*l],sN[2*l+1]);
            }

            for(int l=0; l<nsig; ++l)
            for(int m=0; m<ndir; ++m)
            {
            const int b = g.bin(l,m);
            const double cs = g.costh[m], sn = g.sinth[m];
            const double ns = double(Ns[b]), ms = double(Ms[b]);
            // the face coefficient: mean Sigma along the axis times the energy-weighted age of the two cells
            auto face = [&](const nb &o, const std::vector<double> &so, bool xaxis) -> double
            {
            if(o.N==nullptr)
            return 0.0;
            const double ss = 0.5*(sC[2*l] + so[2*l]), nn = 0.5*(sC[2*l+1] + so[2*l+1]);
            const double dax = xaxis ? ss*cs*cs + nn*sn*sn : ss*sn*sn + nn*cs*cs;
            const double nsum = ns + double(o.Ns[b]);
            const double ta = (nsum>0.0) ? (ms + double(o.Ms[b]))/nsum : 0.0;
            return dt*dax*ta*(xaxis ? rx2 : ry2);
            };
            const double aW = face(W,sW,true), aE = face(E,sE,true), aS = face(S,sS,false), aN = face(Nn,sN,false);
            const double den = 1.0 + aW + aE + aS + aN;
            double num = ns, numM = ms;
                if(W.N)  {num += aW*double(W.N[b]);  numM += aW*double(W.M[b]);}
                if(E.N)  {num += aE*double(E.N[b]);  numM += aE*double(E.M[b]);}
                if(S.N)  {num += aS*double(S.N[b]);  numM += aS*double(S.M[b]);}
                if(Nn.N) {num += aN*double(Nn.N[b]); numM += aN*double(Nn.M[b]);}
                if(cr && ns>0.0)
                {
                const double tc = ms/ns;
                const double dxy = (sC[2*l] - sC[2*l+1])*sn*cs*tc;
                num  += dt*2.0*dxy*rxy*(double(sNE[b]) - double(sNW[b]) - double(sSE[b]) + double(sSW[b]));
                numM += dt*2.0*dxy*rxy*(double(mNE[b]) - double(mNW[b]) - double(mSE[b]) + double(mSW[b]));
                }
            N[b] = float(std::max(num,0.0)/den);
            M[b] = float(std::max(numM,0.0)/den);
            }
        }
    }

    pex->start(p,pgc,*e->N);
    pex->start(p,pgc,*gM);
}
