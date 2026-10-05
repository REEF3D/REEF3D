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
#include"seastate_vtp.h"
#include"seastate_exchange.h"
#include"seastate_implicit.h"
#include"seastate_source.h"
#include"regression_dump.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include<iostream>
#include<iomanip>

seastate_f::seastate_f(lexer *p, ghostcell *pgc) : pprint(nullptr), preg(nullptr), pex(nullptr), psolv(nullptr), N0(nullptr), psrc(nullptr),
                                                  iter_max(1), iter_done(0), dtw(0.0), coupled(false), conv(0.0),
                                                  etot(0.0), hsmax(0.0), hsmean(0.0), nmin(0.0), cells_active(0.0),
                                                  starttime(0.0), endtime(0.0)
{
    side[0]=side[1]=side[2]=side[3]=0;

    e = new fdm_seastate(p);
}

seastate_f::~seastate_f()
{
    delete pprint;
    delete preg;
    delete pex;
    delete psolv;
    delete N0;
    delete psrc;
    delete e->cg;
    delete e->kw;
    delete e->N;
    delete e->grid;
    delete e;
}

void seastate_f::start(lexer *p, ghostcell *pgc)
{
    ini(p,pgc);

    preg = new regression_dump(p);
    preg->seastate_ini(p,e,pgc);

    if(p->mpirank==0)
    cout<<"starting mainloop.SEASTATE"<<endl;

//-----------MAINLOOP SEASTATE----------------------------
    while(p->count<p->N45 && p->simtime<p->N41)
    {
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"------------------------------------"<<endl;
        cout<<p->count<<endl;
        cout<<"simtime: "<<p->simtime<<endl;
        cout<<"timestep: "<<p->dt<<endl;
        }

        step(p,pgc);

        p->simtime+=p->dt;

        // 2D spectra for the wave generation of FNPF/NHFLOW
        if(p->A760>0)
        handover(p,pgc);

        // printer
        double ptime=pgc->timer();
        pprint->start(p,e,pgc);
        log_step(p);
        preg->seastate_step(p,e,pgc);
        p->printouttime=pgc->timer()-ptime;

        endtime=pgc->timer();
        p->itertime=endtime-starttime;
        p->totaltime+=p->itertime;
        p->meantime=p->totaltime/double(p->count);

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"Hs max: "<<setprecision(5)<<hsmax<<"   Hs mean: "<<setprecision(5)<<hsmean<<"   iterations: "<<iter_done;
        if(p->A700==2)
        cout<<"   max. relative change of Hs: "<<setprecision(3)<<conv;
        cout<<endl;
        cout<<"printouttime: "<<setprecision(3)<<p->printouttime<<endl;
        cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
        }

        p->gctime=0.0;
        p->xtime=0.0;
        p->printouttime=0.0;

        if(hsmax!=hsmax || nmin<0.0)
        {
            if(p->mpirank==0)
            cout<<endl<<"EMERGENCY STOP  --  SEASTATE: wave height is NaN or the action density negative"<<endl<<endl;

            pprint->print2D(p,e,pgc);
            pgc->final(true);
        }
    }

    if(p->mpirank==0)
    {
    cout<<endl<<"******************************"<<endl<<endl;

    if(p->plog)
    p->plog->end(p,"finished");

    cout<<"modelled time: "<<p->simtime<<endl;
    cout<<endl;
    }

    preg->seastate_final(p,e,pgc);

    pgc->final();
}

void seastate_f::step(lexer *p, ghostcell *pgc)
{
    // implicit transport in x, y, sigma, theta with the source terms (if any)

    transport(p,pgc);

    parameters(p,pgc);
}

void seastate_f::step_coupled(lexer *p, ghostcell *pgc, double dt)
{
    dtw = dt;

    // the host has written depth, eta, U, V: cells of the initial active set dry and re-wet
    // with the host's water depth; a cell that dries loses its spectrum
    const int nbin = e->grid->nbin;

    IMALOOP
    JMALOOP
    {
    const int w = (e->wet0(i,j)==1 && e->depth(i,j)>=p->A705) ? 1 : 0;

        if(w==0 && e->wet(i,j)==1)
        {
        float *s = e->N->spec(i,j);

        if(s!=nullptr)
        for(int b=0; b<nbin; ++b)
        s[b] = 0.0f;
        }

    e->wet(i,j) = w;
    }

    kinematics(p,pgc,dt);

    transport(p,pgc);

    parameters(p,pgc);

    log_step(p);

    if(p->A760>0)
    handover(p,pgc);
}
