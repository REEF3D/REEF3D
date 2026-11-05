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


#include"seastate_roller.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_exchange.h"
#include"seastate_source.h"
#include"seastate_dispersion.h"
#include"fdm_seastate.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

namespace
{
    inline double pos(double a) {return a>0.0 ? a : 0.0;}
    inline double neg(double a) {return a<0.0 ? a : 0.0;}
}

seastate_roller::seastate_roller(lexer *p, fdm_seastate *e, double beta_) : R(p), Dr(p), Dw(p), beta(beta_)
{
    const seastate_grid &g = *e->grid;

    ndir = g.ndir;

    Rt = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,ndir,p->A704);
    Rt->build(e->wet0.V);
    Rt->fill(0.0f);

    R0 = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,ndir,p->A704);
    R0->build(e->wet0.V);

    Sw = new seastate_store(p->imin,p->jmin,p->imax,p->jmax,ndir,p->A704);
    Sw->build(e->wet0.V);

    pex = new seastate_exchange(p,ndir,1);

    m0.assign(4,ndir);
    m1.assign(4,-1);

    for(int m=0; m<ndir; ++m)
    {
    m0[g.quad[m]] = std::min(m0[g.quad[m]],m);
    m1[g.quad[m]] = std::max(m1[g.quad[m]],m);
    }

    P.assign(g.nbin,0.0);
    D.assign(g.nbin,0.0);

    IMALOOP
    JMALOOP
    R(i,j)=Dr(i,j)=Dw(i,j)=0.0;
}

seastate_roller::~seastate_roller()
{
    delete pex;
    delete Sw;
    delete R0;
    delete Rt;
}

const float *seastate_roller::spec(int i, int j) const
{
    return Rt->spec(i,j);
}

double seastate_roller::celerity(fdm_seastate *e, int i, int j)
{
    const float *k = e->kw->spec(i,j);

    if(k==nullptr || !(k[0]>0.0f))
    return 0.0;

    return e->grid->sig[0]/double(k[0]);
}

void seastate_roller::step(lexer *p, ghostcell *pgc, fdm_seastate *e, seastate_source *src, double dt, int iterations)
{
    const seastate_grid &g = *e->grid;
    const double grav = seastate_gravity;

    R0->copy_from(*Rt);

    // source: breaking dissipation of the waves per direction, from the spectrum after the wave step
    IMALOOP
    JMALOOP
    {
    Dw(i,j) = 0.0;

    float *s = Sw->spec(i,j);

        if(s==nullptr)
        continue;

        for(int m=0; m<ndir; ++m)
        s[m] = 0.0f;

        if(e->wet(i,j)==0 || src==nullptr)
        continue;

    const float *N = e->N->spec(i,j);
    src->compute(N,e->depth(i,j),e->kw->spec(i,j),e->cg->spec(i,j),P.data(),D.data());

        for(int m=0; m<ndir; ++m)
        {
        const double Ew = g.sig[0]*double(N[g.bin(0,m)])*g.dsig[0];
        s[m] = float(src->brk_rate*Ew);
        Dw(i,j) += src->brk_rate*Ew*g.dtheta;
        }
    }

    const double rdt = 1.0/dt;

    for(int it=0; it<std::max(iterations,1); ++it)
    for(int q=0; q<4; ++q)
    {
    sweep(p,e,q,rdt);
    pex->start(p,pgc,*Rt);
    }

    // integrals
    IMALOOP
    JMALOOP
    {
    R(i,j) = Dr(i,j) = 0.0;

    const float *r = Rt->spec(i,j);

        if(r==nullptr || e->wet(i,j)==0)
        continue;

    const double c = celerity(e,i,j);
    const double kd = c>0.0 ? 2.0*grav*beta/c : 0.0;

        for(int m=0; m<ndir; ++m)
        {
        R(i,j)  += double(r[m])*g.dtheta;
        Dr(i,j) += kd*double(r[m])*g.dtheta;
        }
    }
}

void seastate_roller::sweep(lexer *p, fdm_seastate *e, int q, double rdt)
{
    if(m1[q]<m0[q])
    return;

    const seastate_grid &g = *e->grid;
    const double grav = seastate_gravity;

    const bool idown = (q==1 || q==2);
    const bool jdown = (q==2 || q==3);

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
    i = idown ? p->knox-1-ii : ii;
    j = jdown ? p->knoy-1-jj : jj;

    float *r = Rt->spec(i,j);

        if(r==nullptr)
        continue;

        if(e->wet(i,j)==0)
        {
        for(int m=m0[q]; m<=m1[q]; ++m)
        r[m] = 0.0f;
        continue;
        }

    const int ic=i, jc=j;
    const double c = celerity(e,ic,jc);
    const double U = e->U(ic,jc), V = e->V(ic,jc);
    const double kd = c>0.0 ? 2.0*grav*beta/c : 0.0;
    const double rdx = 1.0/p->DXN[IP];
    const double rdy = 1.0/p->DYN[JP];

    const float *r0 = R0->spec(ic,jc);
    const float *sw = Sw->spec(ic,jc);

    // neighbours: wave-active cells only (no roller enters through land, dry cells or the sides)
    const bool wW = e->wet(ic-1,jc)==1, wE = e->wet(ic+1,jc)==1, wS = e->wet(ic,jc-1)==1, wN = e->wet(ic,jc+1)==1;
    const double cW = wW ? celerity(e,ic-1,jc) : 0.0, cE = wE ? celerity(e,ic+1,jc) : 0.0;
    const double cS = wS ? celerity(e,ic,jc-1) : 0.0, cN = wN ? celerity(e,ic,jc+1) : 0.0;
    const float *rW = wW ? Rt->spec(ic-1,jc) : nullptr, *rE = wE ? Rt->spec(ic+1,jc) : nullptr;
    const float *rS = wS ? Rt->spec(ic,jc-1) : nullptr, *rN = wN ? Rt->spec(ic,jc+1) : nullptr;

        for(int m=m0[q]; m<=m1[q]; ++m)
        {
        const double cs = g.costh[m], sn = g.sinth[m];

        const double cxc = c*cs + U, cyc = c*sn + V;
        const double cxw = wW ? 0.5*(cxc + cW*cs + e->U(ic-1,jc)) : cxc;
        const double cxe = wE ? 0.5*(cxc + cE*cs + e->U(ic+1,jc)) : cxc;
        const double cys = wS ? 0.5*(cyc + cS*sn + e->V(ic,jc-1)) : cyc;
        const double cyn = wN ? 0.5*(cyc + cN*sn + e->V(ic,jc+1)) : cyc;

        const double dg = rdt + (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy + kd;

        double rhs = rdt*double(r0[m]) + double(sw[m]);

        if(rW) rhs += pos(cxw)*rdx*double(rW[m]);
        if(rE) rhs -= neg(cxe)*rdx*double(rE[m]);
        if(rS) rhs += pos(cys)*rdy*double(rS[m]);
        if(rN) rhs -= neg(cyn)*rdy*double(rN[m]);

        r[m] = float(rhs/dg);
        }

    i = ic;
    j = jc;
    }
}
