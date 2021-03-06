/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"driver.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"6DOF_header.h"
#include"lexer.h"


void driver::loop_cfd_df(fdm* a)
{
    momentum_RK2_df* pmom_df2 = new momentum_RK2_df(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow);      
    momentum_RK3_df* pmom_df = new momentum_RK3_df(p,a,pgc,pconvec,pdiff,ppress,ppois,pturb,psolv,ppoissonsolv,pflow); 

    sixdof_df* p6dof_df = new sixdof_df(p,a,pgc);
    p6dof_df->initialize(p, a, pgc, pnet);

	driver_ini();

    if(p->mpirank==0)
    cout<<"starting mainloop.CFD 6DOF"<<endl;


//-----------MAINLOOP CFD FSI----------------------------
	while(p->count<p->N45 && p->simtime<p->N41  && p->sedtime<p->S19)
	{
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
            cout<<"------------------------------"<<endl;
            cout<<p->count<<endl;

            cout<<"simtime: "<<p->simtime<<endl;
            cout<<"timestep: "<<p->dt<<endl;

            if(p->B90>0 && p->B92<=11)
            cout<<"t/T: "<<p->simtime/p->wT<<endl;
            
            if(p->B90>0 && p->B92>11)
            cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }

        pflow->flowfile(p,a,pgc,pturb);

        pflow->wavegen_precalc(p,pgc);

        fill_vel(p,a,pgc);
        
        // Benchmark cases
        pbench->start(p,a,pgc,pconvec);
        
        // Free-surface computation
        pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,ppart,a->phi);
        poneph->update(p,a,pgc,pflow);
        
        // Turbulence computation
        pturb->start(a,p,pturbdisc,pturbdiff,psolv,pgc,pflow,pvrans);
        
        // Heat computation
        pheat->start(a,p,pconvec,pdiff,psolv,pgc,pflow);
        
        // Concentration computation
        pconc->start(a,p,pconcdisc,pconcdiff,pturb,psolv,pgc,pflow);
        
        // Suspension computation
        psusp->start(a,p,pconcdisc,psuspdiff,psolv,pgc,pflow);
        
        // Sediment computation
        psed->start(p,a,pconvec,pgc,pflow,ptopo,preto,psusp,pbed);

        pflow->u_relax(p,a,pgc,a->u);
		pflow->v_relax(p,a,pgc,a->v);
		pflow->w_relax(p,a,pgc,a->w);
		pfsf->update(p,a,pgc,a->phi);
	
        if (p->Y2 == 1)
        {
            // Momentum and 6DOF motion
            pmom_df2->starti(p,a,pgc,p6dof_df,pvrans,pnet);

            // Save previous timestep
            pmom_df2->utimesave(p,a,pgc);
            pmom_df2->vtimesave(p,a,pgc);
            pmom_df2->wtimesave(p,a,pgc);
            pflow->veltimesave(p,a,pgc,pvrans);
            pturb->ktimesave(p,a,pgc);
            pturb->etimesave(p,a,pgc);
        }
        else
        {
            // Momentum and 6DOF motion
            pmom_df->starti(p,a,pgc,p6dof_df,pvrans,pnet);

            // Save previous timestep
            pmom_df->utimesave(p,a,pgc);
            pmom_df->vtimesave(p,a,pgc);
            pmom_df->wtimesave(p,a,pgc);
            pflow->veltimesave(p,a,pgc,pvrans);
            pturb->ktimesave(p,a,pgc);
            pturb->etimesave(p,a,pgc);
        }

        //timestep control
        ptstep->start(a,p,pgc,pturb);
        p->simtime+=p->dt;


        // printer
        pprint->start(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,psed);

        // Shell-Printout
        if(p->mpirank==0)
        {
        endtime=pgc->timer();

		p->itertime=endtime-starttime;
		p->totaltime+=p->itertime;
		p->gctotaltime+=p->gctime;
		p->Xtotaltime+=p->xtime;
		p->meantime=(p->totaltime/double(p->count));
		p->gcmeantime=(p->gctotaltime/double(p->count));
		p->Xmeantime=(p->Xtotaltime/double(p->count));

            if( (p->count%p->P12==0))
            {
            if(p->B90>0)
            cout<<"wavegentime: "<<setprecision(3)<<p->wavetime<<endl;

            cout<<"reinitime: "<<setprecision(3)<<p->reinitime<<endl;
            cout<<"gctime: "<<setprecision(3)<<p->gctime<<"\t average gctime: "<<setprecision(3)<<p->gcmeantime<<endl;
            cout<<"Xtime: "<<setprecision(3)<<p->xtime<<"\t average Xtime: "<<setprecision(3)<<p->Xmeantime<<endl;
            cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
            cout<<"timer per step: "<<setprecision(3)<<p->itertime<<endl;
            }
        // Write log files
        mainlog(p);
        maxlog(p);
        solverlog(p);
        if(p->count%p->S44==0 && p->count>=p->S43 && p->S10>0)
        sedimentlog(p);
        }
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;

    stop(p,a,pgc);
	}

	if(p->mpirank==0)
	{
	cout<<endl<<"******************************"<<endl<<endl;

	cout<<"modelled time: "<<p->simtime<<endl;
	cout << endl;

    mainlogout.close();
    maxlogout.close();
    solvlogout.close();
    sedlogout.close();
	}

    pgc->final();
}

