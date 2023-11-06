/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
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
#include"FSI_header.h"

void driver::loop_cfd_df(fdm* a)
{
    if(p->mpirank==0)
    cout<<"starting mainloop.CFD_DF"<<endl;

//-----------MAINLOOP CFD FSI----------------------------
	while(p->count<p->N45 && p->simtime<p->N41  && p->sedtime<p->S19)
	{
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
            cout<<"------------------------------------"<<endl;
            cout<<p->count<<endl;

            cout<<"simtime: "<<setprecision(3)<<p->simtime<<endl;
            cout<<"timestep: "<<p->dt<<endl;
            cout<<"fbtimestep: "<<p->fbdt<<" fbmax: "<<p->fbmax<<endl;
            

            if(p->B90>0 && p->B92<=11)
            cout<<"t/T: "<<p->simtime/p->wT<<endl;
            
            if(p->B90>0 && p->B92>11)
            cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }
        
        p->fbmax=0.0;

        pflow->flowfile(p,a,pgc,pturb);

        pflow->wavegen_precalc(p,pgc);

        fill_vel(p,a,pgc);
        
        // Benchmark cases
        pbench->start(p,a,pgc,pconvec);

        pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,ppls,a->phi);
        poneph->update(p,a,pgc,pflow);

        // Turbulence computation
        pturb->start(a,p,pturbdisc,pturbdiff,psolv,pgc,pflow,pvrans);
        
        // Heat computation
        pheat->start(a,p,pconvec,pdiff,psolv,pgc,pflow);
        
        // Concentration computation
        pconc->start(a,p,pconcdisc,pconcdiff,pturb,psolv,pgc,pflow);
        
        // Sediment computation
        psed->start_cfd(p,a,pgc,pflow,preto,psolv);

        pflow->u_relax(p,a,pgc,a->u);
        pflow->v_relax(p,a,pgc,a->v);
        pflow->w_relax(p,a,pgc,a->w);
        pfsf->update(p,a,pgc,a->phi);
	
        // Momentum and 6DOF motion
        pmom_df->starti(p,a,pgc,p6dof_df,pvrans,pnet,pfsi);

        // Save previous timestep
        pmom_df->utimesave(p,a,pgc);
        pmom_df->vtimesave(p,a,pgc);
        pmom_df->wtimesave(p,a,pgc);
        pflow->veltimesave(p,a,pgc,pvrans);
        pturb->ktimesave(p,a,pgc);
        pturb->etimesave(p,a,pgc);

        //timestep control
        p->simtime+=p->dt;
        ptstep->start(a,p,pgc,pturb);

        // printer
        pprint->start(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);

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
            cout<<"fbtime: "<<setprecision(3)<<p->fbtime<<endl;
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
        }
    p->utime=p->vtime=p->wtime=0.0;
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;
    p->fsitime=0.0;
    p->fbtime=0.0;

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
	}

    pgc->final();
}

