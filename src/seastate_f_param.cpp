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
#include"seastate_amr.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"lexer.h"
#include"ghostcell.h"
#include"seastate_implicit.h"
#include<algorithm>
#include<functional>
#include<vector>

void seastate_f::parameters(lexer *p, ghostcell *pgc)
{
    // integrated wave parameters of every active cell; the rows shared out to the threads (A 798, Phase 13), the totals
    // summed afterwards in the order of the cells, as before
    const int nbin = e->grid->nbin;
    const int ni = p->knox+2*p->margin, nj = p->knoy+2*p->margin;
    const int nt = std::max(1,std::min(p->A798,ni));
    std::vector<double> em(size_t(ni)*nj,0.0), hc(size_t(ni)*nj,0.0);
    std::vector<char> in(size_t(ni)*nj,0);
    std::vector<double> hmt(nt,0.0), vmt(nt,0.0);

    const std::function<void(int)> job = [&](int t)
    {
        seastate_param sp;
        double hm=0.0, vm=0.0;

        for(int i=-p->margin+t; i<p->knox+p->margin; i+=nt)
        for(int j=-p->margin; j<p->knoy+p->margin; ++j)
        {
        e->Hs(i,j)=e->Tm01(i,j)=e->Tm10(i,j)=e->Tp(i,j)=e->dir(i,j)=e->spread(i,j)=0.0;

            if(e->wet(i,j)==1)
            {
            sp.compute(*e->grid,e->N->spec(i,j));

            e->Hs(i,j)     = sp.Hs;
            e->Tm01(i,j)   = sp.Tm01;
            e->Tm10(i,j)   = sp.Tm10;
            e->Tp(i,j)     = sp.Tp;
            e->dir(i,j)    = sp.dir;
            e->spread(i,j) = sp.spread;

                // interior cells: totals
                if(i>=0 && i<p->knox && j>=0 && j<p->knoy && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin)]>0)
                {
                const size_t c = size_t(i+p->margin)*nj + (j+p->margin);
                em[c] = sp.m0*p->DXN[i+increment::marge]*p->DYN[j+increment::marge];
                hc[c] = sp.Hs;
                in[c] = 1;
                hm = std::max(hm,sp.Hs);

                const float *s = e->N->spec(i,j);
                for(int b=0; b<nbin; ++b)
                vm = std::min(vm,double(s[b]));
                }
            }
        }
    hmt[t] = hm;
    vmt[t] = vm;
    };
    seastate_implicit::run_threads(nt,job);

    double esum=0.0, hssum=0.0;
    for(size_t c=0; c<in.size(); ++c)
    if(in[c])
    {
    esum  += em[c];
    hssum += hc[c];
    }
    double hmax = *std::max_element(hmt.begin(),hmt.end());
    double vmin = *std::min_element(vmt.begin(),vmt.end());

    // mesh refinement: parameters of the patches; Hs max and N min over the leaf cells of all grids
    // (E_tot and the mean Hs from level 0, whose covered cells hold the restricted spectra)
    if(pamr!=nullptr)
    pamr->parameters(nullptr,hmax,vmin);

    etot   = pgc->globalsum(esum);
    hsmax  = pgc->globalmax(hmax);
    nmin   = pgc->globalmin(vmin);
    hsmean = cells_active>0.0 ? pgc->globalsum(hssum)/cells_active : 0.0;
}
