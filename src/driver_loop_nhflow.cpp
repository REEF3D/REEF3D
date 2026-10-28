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

#include"driver.h"
#include"dem.h"
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
#include"nhflow_header.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"seastate_nhflow.h"
#include"nhflow_amr.h"
#include"nhflow_longrun.h"
#include"regression_dump.h"
#include"coral/reef3d/coral_reef3d.h"
#include"runlog.h"
#include"lagoon_output.h"
#include"reefrun_ckpt.h"
#include"reefcouple.h"
#include"wave_lib_fnpf_couple.h"

void driver::loop_nhflow()
{
    if(p->mpirank==0)
    cout<<"starting mainloop.NHFLOW"<<endl;
    
//-----------MAINLOOP NHFLOW----------------------------
	while(p->count<p->N45 && p->simtime<p->N41  && p->sedtime<p->S19)
	{		
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"------------------------------------"<<endl;
        cout<<p->count<<endl;
        
        cout<<"simtime: "<<setprecision(6)<<p->simtime<<endl;
        cout<<"timestep: "<<setprecision(4)<<p->dt<<endl;
        
		if(p->B90>0 && p->B92<=11)
		cout<<"t/T: "<<p->simtime/p->wT<<endl;
        
        if(p->B90>0 && p->B92>11)
		cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }
        
        pflow->flowfile(p,a,pgc,pturb);
        
        // runtime coupling (B 92 62): feedback of the bodies to FNPF at the start of the step (I 243)
        if(reefcouple::recv!=nullptr)
        reefcouple::recv->feedback(p,pgc,p6dof,nullptr);
        
        pflow->wavegen_precalc_nhflow(p,d,pgc);
        
        // REEF3D::SEASTATE: wave step and wave forcing (A 750 1)
        if(d->wave!=nullptr)
        d->wave->start(p,d,pgc);
			
        pnhfturb->start(p,d,pgc,pnhfscalarconvec,pnhfturbdiff,psolv,pflow,pnhfvrans);        
        
		// Sediment Computation
        psed->start_susp_nhflow(p,d,pgc,pflow,psolv);
        psed->start_nhflow(p,d,pgc,pflow);
        pnhfsf->depth_update(p,d,pgc,pflow);
        
        // Lagrangian particles: release, velocity at t^n
        pnhfpart->step_begin(p,d,pgc);

        // DEM
        pdem->start_nhflow(p,d,pgc);
        
        // FEM strong coupling (Z 30): the momentum step is repeated until the motion of the
        // rigid bodies converges (no repetition without it)
        pnhfdf->fem_step_begin(p,d,pgc);
        do
        {
        pnhfmom->start(p,d,pgc,pflow,pss,precon,pnhfconvec,pnhfdiff,
                       pnhpress,ppoissonsolv,psolv,pnhf,pnhfsf,pnhfturb,pnhfvrans); 
        }
        while(pnhfdf->fem_step_redo(p,d,pgc));

        //save previous timestep
        pnhfturb->ktimesave(p,d,pgc);
        pnhfturb->etimesave(p,d,pgc);
        //pflow->veltimesave(p,a,pgc,pvrans);
        
        // REEF3D::CORAL (Z 40): structural model from simtime to simtime + dt
        if(pcoral!=nullptr)
        pcoral->step(p,pgc);

        //timestep control
        p->simtime+=p->dt;

        // long-run statistics (P 160), weight dt of the step just done
        if(plongrun!=nullptr)
        plongrun->step(p,d,pgc,pnhfamr);
        
        // Lagrangian particles: Heun corrector with the new field, exchange, output
        pnhfpart->step_end(p,d,pgc);
        
        pnhfstep->start(p,d,pgc);
        
        // printer
        pprint->start(p,d,pgc,pflow,pnhfturb,psed);
        lagoon_output::poll_all(p,pgc);  // outputs of the LAGOON store written by now (P 18)
        preg->nhflow_step(p,d,pgc);
        
        if(pnhfamr!=nullptr)
        pnhfamr->print(p,d,pgc);

        // Shell-Printout
        p->wavecalctime = pgc->globalmax(p->wavecalctime);
        
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
		
		
        if(p->count%p->P12==0)
        {
        if(p->B90>0)
		cout<<"wavegentime: "<<setprecision(5)<<p->wavecalctime<<endl;
		if(p->X10>0)
        cout<<"fbtime: "<<setprecision(3)<<p->fbtime<<endl;
        cout<<"poissontime: "<<setprecision(3)<<p->poissontime<<"   matrixtime: "<<setprecision(3)<<p->matrixtime<<"   poisson cost: "<<100.0*p->poissontime/p->itertime<<" %"<<endl;
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
        
    volumelog(p);
    
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavecalctime=0.0;
	p->field4time=0.0;
    p->fbtime=0.0;
    p->matrixtime=0.0;
    p->poissontime=0.0;

    // checkpoint on any decomposition (P 48), rebalance on the measured load (G 46): after the logs of the
    // step, so that a restart cuts them back behind this step (reefckpt_open)
    ckpt().step();
	
    stop(p,a,pgc);
	}

	// a checkpoint at the end of the run (P 48): it can be continued with a later end (I 45 1)
	ckpt().finish();

	// the last outputs of the LAGOON store are counted (P 18)
	lagoon_output::finish_all(p,pgc);

	if(p->mpirank==0)
	{
	cout<<endl<<"******************************"<<endl<<endl;

	if(p->plog)
	p->plog->end(p,"finished");

	cout<<"modelled time: "<<p->simtime<<endl;
	cout<<endl;

    mainlogout.close();
    maxlogout.close();
    vollogout.close();
    solvlogout.close();
	}

    preg->nhflow_final(p,d,pgc);

    pgc->final();
}

