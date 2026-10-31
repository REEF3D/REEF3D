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

#include"fnpf_print_scour.h"
#include"scour_gauge.h"
#include"fnpf_sigma_col.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

fnpf_print_scour::fnpf_print_scour(lexer *p, fdm_fnpf *c)
{
    if(p->P83>0)
    psf = new scour_gauge(p,0,"FNPF");

    if(p->P84>0)
    psb = new scour_gauge(p,1,"FNPF");
}

fnpf_print_scour::~fnpf_print_scour()
{
    delete psf;
    delete psb;
}

void fnpf_print_scour::start(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    scour_gauge *sg[2] = {psf,psb};

    for(scour_gauge *g : sg)
    if(g!=nullptr)
    {
        sample(p,c,pgc,g);
        g->step(p,smp.data());
    }
}

void fnpf_print_scour::sample(lexer *p, fdm_fnpf *c, ghostcell *pgc, scour_gauge *g)
{
    const int nv = scour_gauge::nsmp;
    smp.assign(g->num*nv,-1.0e20);

    const int sI = p->jmax*p->kmaxF;
    const int sJ = p->kmaxF;
    const double *const Fi = c->Fi;

    for(n=0;n<g->num;++n)
    {
        i = p->posc_i(g->x[n]);
        j = p->j_dir==1 ? p->posc_j(g->y[n]) : 0;

        // the rank whose cell holds the gauge
        if(i<0 || i>=p->knox || j<0 || j>=p->knoy)
        continue;

        double ub=0.0, vb=0.0;
        const double h = c->WL(i,j);

        k=0;
        const int q = FIJK;

        // the bed node: as the interior nodes of velcalc_sig (4th-order horizontal, central in sigma
        // with the ghost node below the bed)
        if(p->flag7[q]>0 && p->wet[IJ]>0)
        {
            fnpf_sigma_col col;

            if(p->sigx==nullptr)
            col.set(p,c,i,j);

            const double dx4 = (-p->XP[IP2] + 8.0*p->XP[IP1] - 8.0*p->XP[IM1] + p->XP[IM2]);
            const double sgx = p->sigx!=nullptr ? p->sigx[q] : col.sx(p->ZN[KP]);

            ub = (-Fi[q+2*sI] + 8.0*Fi[q+sI] - 8.0*Fi[q-sI] + Fi[q-2*sI])/dx4
               + sgx*((Fi[q+1]-Fi[q-1])/(p->DZN[KP]+p->DZN[KM1]));

            if(p->j_dir==1)
            {
            const double dy4 = (-p->YP[JP2] + 8.0*p->YP[JP1] - 8.0*p->YP[JM1] + p->YP[JM2]);
            const double sgy = p->sigy!=nullptr ? p->sigy[q] : col.sy(p->ZN[KP]);

            vb = (-Fi[q+2*sJ] + 8.0*Fi[q+sJ] - 8.0*Fi[q-sJ] + Fi[q-2*sJ])/dy4
               + sgy*((Fi[q+1]-Fi[q-1])/(p->DZN[KP]+p->DZN[KM1]));
            }
        }

        smp[n*nv+0] = ub;
        smp[n*nv+1] = vb;
        smp[n*nv+2] = ub;   // potential flow: the current is depth uniform, its mean is taken in scour_gauge
        smp[n*nv+3] = vb;
        smp[n*nv+4] = p->wet[IJ]>0 ? h : 0.0;
    }

    pgc->globalmax(smp.data(),g->num*nv);   // one reduction for all gauges
}

void fnpf_print_scour::ckpt_state(std::vector<double> &v) const
{
    v.push_back(psf!=nullptr ? 1.0 : 0.0);
    if(psf!=nullptr)
    psf->ckpt_state(v);

    v.push_back(psb!=nullptr ? 1.0 : 0.0);
    if(psb!=nullptr)
    psb->ckpt_state(v);
}

void fnpf_print_scour::ckpt_restore(const std::vector<double> &v, size_t &pos)
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
