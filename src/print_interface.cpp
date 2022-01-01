/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"print_interface.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"heat.h"
#include"solver.h"
#include"print_wsf.h"
#include"print_wsf_theory.h"
#include"print_wsfline.h"
#include"print_wsfline_y.h"
#include"vorticity_f.h"
#include"vorticity_void.h"
#include"probe_point.h"
#include"probe_line.h"
#include"bedprobe_point.h"
#include"bedprobe_max.h"
#include"ioflow.h"
#include"data.h"
#include"concentration.h"
#include"gage_discharge.h"
#include"fsf_vtp.h"
#include"state.h"
#include"bedshear_probe.h"
#include"sediment.h"
#include"sloshing_force.h"
#include"print_porous.h"
#include<sys/stat.h>
#include<sys/types.h>

print_interface::print_interface(lexer* p, fdm *a, ghostcell *pgc) : nodefill(p), eta(p)
{
	
    if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;
	
	if(p->F50==1)
	gcval_phiext=61;

	if(p->F50==2)
	gcval_phiext=62;

	if(p->F50==3)
	gcval_phiext=63;

	if(p->F50==4)
	gcval_phiext=64;

    if(p->I40==0)
    {
	p->printtime=0.0;
	p->sedprinttime=0.0;
	p->fsfprinttime=0.0;
	p->probeprinttime=0.0;
	p->stateprinttime=0.0;
	p->sedsimtime = 0.0;
	p->sedwavetime = 0.0;
    }
	
	p->Darray(printtime_wT,p->P35);
	
	for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn];
	
	pwsf=new print_wsf(p,a,pgc,0);
	pwsf_theory=new print_wsf_theory(p,a,pgc,0);
	pwsfline=new print_wsfline(p,a,pgc);
	pwsfline_y=new print_wsfline_y(p,a,pgc);
	pprobe = new probe_point(p,a,pgc);
	pline = new probe_line(p,a,pgc);
	pq = new gage_discharge(p,a,pgc);
	
	if(p->P180==1)
	pfsf = new fsf_vtp(p,a,pgc);
	
	if(p->P75==0)
	pvort = new vorticity_void(p,a);

	if(p->P75==1)
	pvort = new vorticity_f(p,a);

	if(p->P121>0)
	pbedpt = new bedprobe_point(p,a,pgc);
	
	if(p->P122>0)
	pbedmax = new bedprobe_max(p,a,pgc);
	
	if(p->P125>0)
	pbedshear = new bedshear_probe(p,a,pgc);

	if(p->P40>0)
	pstate=new state(p,a,pgc);
    
    if(p->P101>0)
	pslosh=new sloshing_force(p,a,pgc);
	
	if(p->B270>0 || p->B274>0 || p->B281>0 || p->B291>0)
	{
	ppor=new print_porous(p,a,pgc);
	ppor->start(p,a,pgc);
	}

	p->printcount=0;
    
    phase=0.0;
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_VTU",0777);
}

print_interface::~print_interface()
{
}


void print_interface::start(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, sediment *psed)
{		
		
		// Print out based on iteration
        if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)
		{
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,psed);
		}
		
		// Print out based on time
        if(p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1 || (p->count==0 &&  p->P30>0.0))
        {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,psed);
		
        p->printtime+=p->P30;
        }
		
		// Print out based on sediment time
        if(p->simtime>p->sedprinttime && p->P34>0.0 && p->P30<0.0 && p->P10==1 || (p->count==0 &&  p->P34>0.0))
        {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,psed);
		
        p->sedprinttime+=p->P34;
        }
				
		// Print out based on time interval
		if(p->P10==1 && p->P35>0)
		for(int qn=0; qn<p->P35; ++qn)
		if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
		{
		print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,psed);	
			
		printtime_wT[qn]+=p->P35_dt[qn];
		}
		
		
		if((p->P62>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P62>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
        pline->start(p,a,pgc,pturb);
        
		
		if(p->P50>0)
        pwsf_theory->height_gauge(p,a,pgc,pflow,a->phi);
		
        if(p->P51>0)
        pwsf->height_gauge(p,a,pgc,a->phi);
		
		if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
        pwsfline->wsfline(p,a,pgc,pflow);
		
		if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
        pwsfline_y->wsfline(p,a,pgc,pflow);
		
		if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
        {
        p->probeprinttime+=p->P55;
        }

		if(p->P61>0)
        pprobe->start(p,a,pgc,pturb);
		
		if(p->P67>0)
		pq->start(p,a,pgc);
        
        if(p->P101>0)
        pslosh->start(p,a,pgc);
		
		// sediment probes
		if(((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->count>=p->S43) || (p->S41==3 && p->count>=p->S43) ) && p->S10>0)
		if((p->S42==1 && p->count%p->S44==0) || (p->S42==2 && p->simtime>=p->sedsimtime) || (p->S42==3  && p->simtime/p->wT>=p->sedwavetime))
		{
		if(p->P121>0)
        pbedpt->bed_gauge(p,a,pgc);
		
		if(p->P122>0)
        pbedmax->bed_max(p,a,pgc); 
		
		if(p->P125>0)
        pbedshear->bedshear_gauge(p,a,pgc,psed);
		}
    
		
		if(p->count%p->P181==0 && p->P182<0.0 && p->P180==1)
		pfsf->start(p,a,pgc);
		
		
		if((p->simtime>p->fsfprinttime && p->P182>0.0 && p->P180==1) || (p->count==0 &&  p->P182>0.0))
        {
        pfsf->start(p,a,pgc);
        p->fsfprinttime+=p->P182;
        }
		
		// Print state out based on iteration
        if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && p->P41>0)
		{
        pstate->write(p,a,pgc,pturb);
		}
		
		// Print sate out based on time
        if(p->simtime>p->stateprinttime && p->P42>0.0 && p->P40>0 || (p->count==0 &&  p->P42>0.0))
        {
        pstate->write(p,a,pgc,pturb);
		
        p->stateprinttime+=p->P42;
        }
}

void print_interface::print3D(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, sediment *psed)
{
	
	
}
