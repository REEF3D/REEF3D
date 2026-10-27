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
#include"seastate_forcing.h"
#include"met_provider.h"
#include"seastate_param.h"
#include"slice4.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<string>

/*--------------------------------------------------------------------
External forcing (A 711 3, A 730 2), see seastate_f.h and
seastate_forcing.h. The wind field (A 730 2) is that of the REEF3D met
forcing W 200 (met_provider: uniform W 201, met-forcing.r3m, the vortex
met-vortex.trk, the blend W 220-227), the same as the surge of a coupled
SFLOW/NHFLOW run; the 10 m wind at the cell centres, not ramped (W 204
ramps the surface stress of the hydrodynamics only).
--------------------------------------------------------------------*/

namespace
{
    void stop(lexer *p, ghostcell *pgc, const std::string &msg)
    {
        if(p->mpirank==0)
        cout<<endl<<"SEASTATE input error  --  "<<msg<<endl<<endl;

        pgc->final(true);
    }

    std::string date(double t)
    {
        // seconds since 1970 -> YYYY-MM-DD HH:MM:SS (H. Hinnant, civil_from_days)
        long s = long(std::floor(t+0.5));
        long z = (s>=0 ? s : s-86399)/86400;
        const long sod = s - z*86400;
        z += 719468;
        const long era = (z>=0 ? z : z-146096)/146097;
        const long doe = z - era*146097;
        const long yoe = (doe - doe/1460 + doe/36524 - doe/146096)/365;
        const long doy = doe - (365*yoe + yoe/4 - yoe/100);
        const long mp = (5*doy + 2)/153;
        const long d = doy - (153*mp+2)/5 + 1;
        const long m = mp<10 ? mp+3 : mp-9;
        const long y = yoe + era*400 + (m<=2);
        char buf[64];
        snprintf(buf,sizeof(buf),"%04ld-%02ld-%02ld %02ld:%02ld:%02ld",y,m,d,sod/3600,(sod/60)%60,sod%60);
        return buf;
    }
}

void seastate_f::forcing_ini(lexer *p, ghostcell *pgc)
{
    std::string err;
    double tfirst = -1.0;

    if(p->A711==3)
    {
    bser = new seastate_spc_series;

        if(!bser->open("seastate-boundary.spc",*e->grid,err))
        stop(p,pgc,"A 711 3: "+err);

    if(!bser->stationary())
    tfirst = bser->time(0);
    }

    if(p->A730==2)
    {
    wmet = met_provider::get(p,pgc);

    // the box of the cell centres incl. the ghost cells: the part of the met file the rank reads (the patches of
    // the mesh refinement lie inside)
    double xs = 1.0e300, xe = -1.0e300, ys = 1.0e300, ye = -1.0e300;

        IMALOOP
        {
        xs = std::min(xs,p->XP[IP]);
        xe = std::max(xe,p->XP[IP]);
        }

        JMALOOP
        {
        ys = std::min(ys,p->YP[JP]);
        ye = std::max(ye,p->YP[JP]);
        }

    wmet->window(xs,xe,ys,ye);

    wU10 = new slice4(p);
    wdir = new slice4(p);
    wUx  = new slice4(p);
    wUy  = new slice4(p);
    }

    // model time 0: A 780, the first time of the files; with the met forcing that of W 200 (W 203), checked against A 780
    tref = (p->A780>0.0) ? seastate_datetime(p->A780) : std::max(tfirst,0.0);

    // (a uniform wind W 200 1 without W 203 has no date: the files keep theirs)
    if(wmet!=nullptr && (p->W200>=2 || p->W203>0.0))
    {
        if(p->A780>0.0 && std::fabs(seastate_datetime(p->A780) - wmet->tref)>0.5)
        stop(p,pgc,"A 780 and the met forcing W 200: model time 0 differs ("+date(seastate_datetime(p->A780))+" against "+date(wmet->tref)
                   +"); set W 203 to the date-time of A 780 or leave A 780 at 0");

    tref = wmet->tref;
    }

    if(p->mpirank==0)
    {
    cout<<"SEASTATE forcing: model time 0 = "<<date(tref)<<((wmet!=nullptr && (p->W200>=2 || p->W203>0.0)) ? " (the met forcing W 200)" : (p->A780>0.0 ? " (A 780)" : " (first time in the files)"))<<endl;

        if(bser!=nullptr)
        {
        cout<<"  seastate-boundary.spc: "<<bser->nloc<<" locations, "<<bser->f.size()<<" frequencies "<<bser->f.front()<<" - "<<bser->f.back()
            <<" Hz, "<<bser->dir.size()<<" directions, ";
        if(bser->stationary())
        cout<<"stationary"<<endl;
        else
        cout<<"first record "<<date(bser->time(0))<<endl;
        }

        if(wmet!=nullptr)
        cout<<"  wind: the met forcing W 200 "<<p->W200<<" (A 730 2), not ramped"<<endl;
    }
}

void seastate_f::boundary_series(lexer *p, ghostcell *pgc)
{
    const int nbin = e->grid->nbin;

    // boundary cells of this rank: face midpoints and the two nearest locations
    for(int s=0; s<4; ++s)
    {
    Nside[s].clear();
    bw[s].clear();

    const bool mine = (s==0 && p->origin_i==0) || (s==1 && p->origin_i+p->knox==p->gknox)
                   || (s==2 && p->origin_j==0) || (s==3 && p->origin_j+p->knoy==p->gknoy);

        if(side[s]!=1 || !mine)
        continue;

    const int n = (s<2) ? p->knoy : p->knox;
    Nside[s].assign(size_t(n)*nbin,0.0f);
    bw[s].resize(n);

        for(int k=0; k<n; ++k)
        {
        double xb, yb;

            if(s<2)
            {
            j = k;
            i = (s==0) ? 0 : p->knox;
            xb = p->XN[IP];
            yb = p->YP[JP];
            }
            else
            {
            i = k;
            j = (s==2) ? 0 : p->knoy;
            xb = p->XP[IP];
            yb = p->YN[JP];
            }

        int a=0, b=0;
        double da=1.0e300, db=1.0e300;

            for(int l=0; l<bser->nloc; ++l)
            {
            const double d = std::hypot(bser->xs[l]-xb,bser->ys[l]-yb);

                if(d<da)
                {
                db = da; b = a;
                da = d;  a = l;
                }
                else if(d<db)
                {
                db = d;  b = l;
                }
            }

        bweight w;
        w.a = a;
        w.b = (bser->nloc>1) ? b : a;

            if(bser->nloc==1 || da<1.0e-9*(1.0+db))
            {
            w.wa = 1.0f;
            w.wb = 0.0f;
            }
            else
            {
            // inverse distance: linear between two locations on a straight side
            w.wa = float(db/(da+db));
            w.wb = float(da/(da+db));
            }

        bw[s][k] = w;
        }
    }

    if(p->mpirank==0)
    cout<<"SEASTATE boundary (A 711 3): spectra of the boundary cells from the two nearest of "<<bser->nloc<<" locations, linear in time"<<endl<<endl;
}

void seastate_f::forcing_update(lexer *p, ghostcell *pgc, double t)
{
    std::string err;
    const double tf = tref + t;

    if(bser!=nullptr)
    {
        if(!bser->advance(tf,err))
        stop(p,pgc,"A 711 3: "+err);

    const double w = bser->weight(tf);
    const int nbin = e->grid->nbin;

        for(int s=0; s<4; ++s)
        for(size_t k=0; k<bw[s].size(); ++k)
        {
        const bweight &q = bw[s][k];
        const float *a0 = bser->N(0,q.a).data(), *b0 = bser->N(0,q.b).data();
        const float *a1 = bser->N(1,q.a).data(), *b1 = bser->N(1,q.b).data();
        float *out = &Nside[s][k*nbin];

            for(int bb=0; bb<nbin; ++bb)
            out[bb] = float((1.0-w)*(q.wa*a0[bb] + q.wb*b0[bb]) + w*(q.wa*a1[bb] + q.wb*b1[bb]));
        }
    }

    if(wmet!=nullptr)
    {
    // model time t of the met forcing (its own model time 0 = tref of SEASTATE)
    wmet->update(p,pgc,t);

        IMALOOP
        JMALOOP
        {
        double u, v, pa;
        wmet->at(p->XP[IP],p->YP[JP],u,v,pa);

        (*wUx)(i,j)  = u;
        (*wUy)(i,j)  = v;
        (*wU10)(i,j) = std::sqrt(u*u + v*v);
        (*wdir)(i,j) = std::atan2(v,u);
        }

        // mesh refinement: the wind of the patch cells
        if(pamr!=nullptr)
        pamr->wind(tf);
    }
}
