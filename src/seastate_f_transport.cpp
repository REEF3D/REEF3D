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
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_param.h"
#include"seastate_exchange.h"
#include"seastate_implicit.h"
#include"seastate_swan_spc.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

void seastate_f::boundary(lexer *p, ghostcell *pgc)
{
    side[0] = p->A712_xm;
    side[1] = p->A712_xp;
    side[2] = p->A712_ym;
    side[3] = p->A712_yp;

    Nb.clear();

    if(p->A711==1)
    parametric_spectrum(p,pgc,Nb,"A 711 1");

    // SWAN 2D spectrum file, first location and time
    if(p->A711==2)
    {
    seastate_swan_spc spc;
    std::string err;

        if(!spc.read("seastate-boundary.spc",err))
        {
            if(p->mpirank==0)
            cout<<endl<<"SEASTATE input error  --  A 711 2: "<<err<<endl<<endl;

            pgc->final(true);
        }

    spc.to_grid(*e->grid,Nb);

    seastate_param sp;
    sp.compute(*e->grid,Nb.data());

        if(p->mpirank==0)
        cout<<"SEASTATE boundary spectrum (A 711 2): seastate-boundary.spc, "<<spc.f.size()<<" frequencies "<<spc.f.front()<<" - "<<spc.f.back()
            <<" Hz, "<<spc.dir.size()<<" directions; on the spectral grid: Hs "<<sp.Hs<<" m, Tp "<<sp.Tp<<" s, direction "<<sp.dir<<" deg"<<endl<<endl;
    }
}

void seastate_f::transport(lexer *p, ghostcell *pgc)
{
    const bool stationary = (p->A700==2);
    const double rdt = stationary ? 0.0 : 1.0/p->dt;
    const bool refraction = (p->A713==1);
    const bool fshift = (p->A714==1);

    iter_done = 0;
    conv = 0.0;

    if(!stationary)
    {
        if(N0!=nullptr)
        N0->copy_from(*e->N);

        for(int it=0; it<iter_max; ++it)
        {
        psolv->iterate(p,pgc,e,pex,N0,rdt,Nb,side,refraction,fshift);
        ++iter_done;
        }

        return;
    }

    // stationary: iterate until the largest relative change of Hs is below A 708
    const seastate_grid &g = *e->grid;
    seastate_param sp;

    hs_old.assign(size_t(p->knox)*p->knoy,0.0);

    SLICELOOP4
    if(e->wet(i,j)==1)
    {
    sp.compute(g,e->N->spec(i,j));
    hs_old[size_t(i)*p->knoy+j] = sp.Hs;
    }

    for(int it=0; it<iter_max; ++it)
    {
    psolv->iterate(p,pgc,e,pex,nullptr,rdt,Nb,side,refraction,fshift);
    ++iter_done;

    double hmax=0.0;
    std::vector<double> hs(size_t(p->knox)*p->knoy,0.0);

        SLICELOOP4
        if(e->wet(i,j)==1)
        {
        sp.compute(g,e->N->spec(i,j));
        hs[size_t(i)*p->knoy+j] = sp.Hs;
        hmax = std::max(hmax,sp.Hs);
        }

    hmax = pgc->globalmax(hmax);

    double cmax=0.0;
    const double floor = std::max(0.01*hmax,1.0e-6);

        SLICELOOP4
        if(e->wet(i,j)==1)
        {
        const size_t n = size_t(i)*p->knoy+j;
        cmax = std::max(cmax,std::fabs(hs[n]-hs_old[n])/std::max(hs[n],floor));
        }

    conv = pgc->globalmax(cmax);
    hs_old.swap(hs);

        if(conv<p->A708)
        break;
    }
}
