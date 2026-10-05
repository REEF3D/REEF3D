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

#include"spectral_f.h"
#include"fdm_spectral.h"
#include"spectral_grid.h"
#include"spectral_store.h"
#include"spectral_vtp.h"
#include"regression_dump.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include<iostream>
#include<iomanip>

spectral_f::spectral_f(lexer *p, ghostcell *pgc) : pprint(nullptr), preg(nullptr),
                                                  etot(0.0), hsmax(0.0), hsmean(0.0), cells_active(0.0),
                                                  starttime(0.0), endtime(0.0)
{
    e = new fdm_spectral(p);
}

spectral_f::~spectral_f()
{
    delete pprint;
    delete preg;
    delete e->N;
    delete e->grid;
    delete e;
}

void spectral_f::start(lexer *p, ghostcell *pgc)
{
    ini(p,pgc);

    preg = new regression_dump(p);
    preg->spectral_ini(p,e,pgc);

    if(p->mpirank==0)
    cout<<"starting mainloop.Spectral"<<endl;

//-----------MAINLOOP Spectral----------------------------
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

        // printer
        double ptime=pgc->timer();
        pprint->start(p,e,pgc);
        log_step(p);
        preg->spectral_step(p,e,pgc);
        p->printouttime=pgc->timer()-ptime;

        endtime=pgc->timer();
        p->itertime=endtime-starttime;
        p->totaltime+=p->itertime;
        p->meantime=p->totaltime/double(p->count);

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"Hs max: "<<setprecision(5)<<hsmax<<"   Hs mean: "<<setprecision(5)<<hsmean<<endl;
        cout<<"printouttime: "<<setprecision(3)<<p->printouttime<<endl;
        cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
        }

        p->gctime=0.0;
        p->xtime=0.0;
        p->printouttime=0.0;

        if(hsmax!=hsmax)
        {
            if(p->mpirank==0)
            cout<<endl<<"EMERGENCY STOP  --  Spectral: wave height is NaN"<<endl<<endl;

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

    preg->spectral_final(p,e,pgc);

    pgc->final();
}

void spectral_f::step(lexer *p, ghostcell *pgc)
{
    // Phase 1: implicit transport in x, y, sigma, theta
    // Phase 2: source terms

    parameters(p,pgc);
}
