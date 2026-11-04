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
with 8. The coupled form in the sweeps
needs the diffusion converged within the step (validation 79: one
iteration per step lags, the old time level slows the propagation).
--------------------------------------------------------------------*/

void seastate_f::gse_step(lexer *p, ghostcell *pgc)
{
    if(p->A776==0 || p->A700!=1)
    return;

    const double dt = dtw;
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
