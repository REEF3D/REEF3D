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

#include"nhflow_print_scour.h"
#include"scour_gauge.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<algorithm>

nhflow_print_scour::nhflow_print_scour(lexer *p, fdm_nhf *d)
{
    if(p->P83>0)
    psf = new scour_gauge(p,0,"NHFLOW");

    if(p->P84>0)
    psb = new scour_gauge(p,1,"NHFLOW");
}

nhflow_print_scour::~nhflow_print_scour()
{
    delete psf;
    delete psb;
}

void nhflow_print_scour::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    scour_gauge *sg[2] = {psf,psb};

    for(scour_gauge *g : sg)
    if(g!=nullptr)
    {
        sample(p,d,pgc,g);
        g->step(p,smp.data());
    }
}

void nhflow_print_scour::sample(lexer *p, fdm_nhf *d, ghostcell *pgc, scour_gauge *g)
{
    const int nv = scour_gauge::nsmp;
    smp.assign(g->num*nv,-1.0e20);

    for(n=0;n<g->num;++n)
    {
        const double xp = g->x[n];
        double yp = g->y[n];

        i = p->posc_i(xp);
        j = p->j_dir==1 ? p->posc_j(yp) : 0;

        // the rank whose cell holds the gauge
        if(i<0 || i>=p->knox || j<0 || j>=p->knoy)
        continue;

        if(p->j_dir==0)
        yp = p->YP[JP];

        const double h   = p->ccslipol4(d->WL,xp,yp);
        const double zb  = p->ccslipol4(d->bed,xp,yp);

        double ub=0.0, vb=0.0, uc=0.0, vc=0.0;

        if(h>1.0e-6)
        {
            // the bed velocity: the centre of the lowest cell (below it ccipol4V goes towards the bed
            // value, inside the bottom boundary layer), the current: D/2 above the bed, at least there
            k=0;
            const double z1 = p->ZSP[IJK];

            ub = p->ccipol4V(d->U,d->WL,d->bed,xp,yp,z1);
            vb = p->ccipol4V(d->V,d->WL,d->bed,xp,yp,z1);

            const double zc = std::max(zb + std::min(0.5*g->D[n],h), z1);
            uc = p->ccipol4V(d->U,d->WL,d->bed,xp,yp,zc);
            vc = p->ccipol4V(d->V,d->WL,d->bed,xp,yp,zc);
        }

        smp[n*nv+0] = ub;
        smp[n*nv+1] = vb;
        smp[n*nv+2] = uc;
        smp[n*nv+3] = vc;
        smp[n*nv+4] = std::max(h,0.0);
    }

    pgc->globalmax(smp.data(),g->num*nv);   // one reduction for all gauges
}

void nhflow_print_scour::ckpt_state(std::vector<double> &v) const
{
    v.push_back(psf!=nullptr ? 1.0 : 0.0);
    if(psf!=nullptr)
    psf->ckpt_state(v);

    v.push_back(psb!=nullptr ? 1.0 : 0.0);
    if(psb!=nullptr)
    psb->ckpt_state(v);
}

void nhflow_print_scour::ckpt_restore(const std::vector<double> &v, size_t &pos)
{
    scour_gauge *sg[2] = {psf,psb};

    for(scour_gauge *g : sg)
    {
        const bool has = pos<v.size() && v[pos++]>0.5;

        if(has)
        {
            if(g!=nullptr)
            g->ckpt_restore(v,pos);
            else
            {
                // a gauge set of the checkpoint not in this run: skip it
                const int nc = pos<v.size() ? int(v[pos]) : 0;
                pos += 3 + 15*nc;
            }
        }
    }
}
