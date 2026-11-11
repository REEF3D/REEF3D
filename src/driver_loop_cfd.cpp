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
#include"cfd_pocket.h"
#include"regression_dump.h"
#include"coral/reef3d/coral_reef3d.h"
#include"dem.h"
#include"pic.h"
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
#include"runlog.h"
#include"lagoon_output.h"
#include"cfd_amr.h"
#include"reefrun_ckpt.h"
#include"reefcouple.h"
#include"wave_lib_fnpf_couple.h"
#include"soil_stress.h"


namespace
{
// a 6DOF body (X 10 1) with a free degree of freedom (X 11, X 111: 1 free, 2 prescribed, 0 fixed)
bool cfd_free_body(lexer *p)
{
    if(p->X10!=1)
    return false;

    const int nb = MAX(p->X20,1);

    for(int b=0; b<nb; ++b)
    {
        int dof[6] = {p->X11_u,p->X11_v,p->X11_w,p->X11_p,p->X11_q,p->X11_r};

        for(int n=0; n<p->X111; ++n)
        if(p->X111_id[n]==b)
        {
        dof[0] = p->X111_u[n];
        dof[1] = p->X111_v[n];
        dof[2] = p->X111_w[n];
        dof[3] = p->X111_p[n];
        dof[4] = p->X111_q[n];
        dof[5] = p->X111_r[n];
        }

        for(int q=0; q<6; ++q)
        if(dof[q]==1 && !(p->j_dir==0 && (q==1 || q==3 || q==5)))
        return true;
    }

    return false;
}
}

void driver::loop_cfd(fdm* a)
{
    if(p->mpirank==0)
    cout<<"starting mainloop.CFD"<<endl;
    
//-----------MAINLOOP CFD----------------------------
	while(p->count<p->N45 && p->simtime<p->N41  && p->sedtime<p->S19)
	{		
        ++p->count;
        starttime=pgc->timer();
        p->poissontime=0.0;     // the projections of the step add their time (pjm_corr, cfd_amr)
        
        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"------------------------------------"<<endl;
        cout<<p->count<<endl;
        
        cout<<"simtime: "<<p->simtime<<endl;
		cout<<setprecision(5)<<"timestep: "<<p->dt<<endl;
        
        if(p->X10>0)
        cout<<"fbtimestep: "<<p->fbdt<<" fbmax: "<<p->fbmax<<endl;
        
		if(p->B90>0 && p->B92<=11)
		cout<<"t/T: "<<p->simtime/p->wT<<endl;
        
        if(p->B90>0 && p->B92>11)
		cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }
        
        pflow->flowfile(p,a,pgc,pturb);
        
        // runtime coupling (B 92 62): feedback to FNPF at the start of the step (I 243, I 244)
        if(reefcouple::recv!=nullptr)
        reefcouple::recv->feedback(p,pgc,p6dof,a);
        
        pflow->wavegen_precalc(p,pgc);

		    
            pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,a->phi);
            pturb->start(a,p,pturbdisc,pturbdiff,psolv,pgc,pflow,pvrans);
            // with the mesh refinement it transports heat and concentration on all grids
            if(pcfdamr==nullptr || !pcfdamr->scalars())
            {
            pheat->start(a,p,pheatdisc,pdiff,psolv,pgc,pflow);
            pconc->start(a,p,pconcdisc,pconcdiff,pturb,psolv,pgc,pflow);
            }
            pmp->start(p,a,pgc,pmpconvec,psolv,pflow,preini);
        
        psed->start_susp(p,a,pgc,pflow,psolv);
        psed->start_cfd(p,a,pgc,pflow,preto,psolv);
        if(pcfdamr!=nullptr)
        pcfdamr->cpm_step(p,a,pgc);
        
        pflow->u_relax(p,a,pgc,a->u);
        pflow->v_relax(p,a,pgc,a->v);
        pflow->w_relax(p,a,pgc,a->w);
        pfsf->update(p,a,pgc,a->phi);
        pdem->start_cfd(p,a,pgc);
        ppic->step_begin(p,a,pgc);   // PIC (F 81 2): the particle velocities to the water faces
        // F 41 2: the pressures of the closed air pockets from their volumes (cfd_pocket.h)
        cfd_pocket_update(p,a,pgc,pmom,pcfdamr);
        if(pcfdamr!=nullptr)
        pcfdamr->step(p,a,pgc,pvrans,p6dof);
        else
        pmom->start(p,a,pgc,pvrans,p6dof); 
        pbench->start(p,a,pgc,pconvec);
		
        //save previous timestep
        pturb->ktimesave(p,a,pgc);
        pturb->etimesave(p,a,pgc);
        pflow->veltimesave(p,a,pgc,pvrans);
        psed->ctimesave(p,a);

        // REEF3D::CORAL (Z 40): structural model from simtime to simtime + dt
        if(pcoral!=nullptr)
        pcoral->step(p,pgc);

        //timestep control
        p->simtime+=p->dt;

        // PIC particles (F 81): moved over the step with the new velocity, exchanged, output
        ppic->start(p,a,pgc);

        // soil stress model (S 200), at the new time level
        if(psoil!=nullptr)
        psoil->step(p,a,pgc);

        ptstep->start(a,p,pgc,pturb);

        // floating bodies with free DOFs and the one-step scheme (N 40 5): the loads of p(n) answer
        // the velocity change of the body over the last step, divided by the new step. A step much
        // longer than the last one (after the first steps, which the forcing limits) turns the loads
        // of the short step into a kick of the body: the step grows by 10 % per step at most
        if(p->N40==5 && p->N48>0 && cfd_free_body(p))
        p->dt = MIN(p->dt, 1.1*p->dt_old);

        if(pcfdamr!=nullptr)
        pcfdamr->timestep(p,a,pgc);
        
        // printer
        if(pcfdamr!=nullptr)
        pcfdamr->print(p,pgc);
        pprint->start(p,a,pgc,pturb,pheat,pflow,pdata,pconc,pmp,psed);
        lagoon_output::poll_all(p,pgc);  // outputs of the LAGOON store written by now (P 18)
        preg->cfd_step(p,a,pgc,pturb,pconc);

        // Shell-Printout
        if(p->mpirank==0)
        {
        endtime=pgc->timer();
        
		p->itertime=endtime-starttime;
		p->totaltime+=p->itertime;
		p->gctotaltime+=p->gctime;
		p->Xtotaltime+=p->xtime;
		p->poissontotaltime+=p->poissontime;
		p->meantime=(p->totaltime/double(p->count));
		p->gcmeantime=(p->gctotaltime/double(p->count));
		p->Xmeantime=(p->Xtotaltime/double(p->count));
        
            if( (p->count%p->P12==0))
            {
            if(p->B90>0)
            cout<<"wavegentime: "<<setprecision(3)<<p->wavecalctime<<endl;
            if(p->X10>0)
            cout<<"fbtime: "<<setprecision(3)<<p->fbtime<<endl;
            cout<<"reinitime: "<<setprecision(3)<<p->reinitime<<endl;
            // all Poisson solves of the step (every stage; with the mesh refinement every projection of all levels)
            cout<<"poissontime: "<<setprecision(3)<<p->poissontime<<"\t total poissontime: "<<setprecision(6)<<p->poissontotaltime
                <<"   poisson cost: "<<setprecision(3)<<(p->itertime>0.0 ? 100.0*p->poissontime/p->itertime : 0.0)<<" %"<<endl;
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
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavecalctime=0.0;
	p->field4time=0.0;
    
    pgc->gcparax(p,a->press,4);

    // checkpoint on any decomposition (P 48), rebalance on the measured load (G 46)
    ckpt().step();
    
    stop(p,a,pgc);
	}

	// the state at the end of the run (N 45, N 41, S 19), unless the last step wrote it: before the
	// checkpoint, which keeps the output counters
	pprint->print_end(p,a,pgc,pturb,pheat,pflow,pdata,pconc,pmp,psed);

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
	cout << endl;

    mainlogout.close();
    maxlogout.close();
    solvlogout.close();
	}

    preg->cfd_final(p,a,pgc,pturb,pconc);

    pgc->final();
}

