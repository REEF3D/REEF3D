/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"sflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"iowave.h"
#include"sediment.h"
#include"hypre_struct2D.h"
#include"sflow_etimestep.h"
#include"sflow_weno_flux.h"
#include"sflow_eta.h"
#include"sflow_momentum.h"
#include"sflow_hydrostatic.h"
#include"sflow_vtp_fsf.h"
#include"sflow_vtp_bed.h"
#include"sflow_filter.h"
#include"sflow_turbulence.h"
#include"6DOF_sflow.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

sflow_f::sflow_f(lexer *p, fdm2D *b, ghostcell* pgc, patchBC_interface *ppBC)
{
    pBC = ppBC;
    
	maxcoor(p,b,pgc);
}

sflow_f::~sflow_f()
{
}

void sflow_f::start(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	
	// logic
	logic(p,b,pgc);
	
	// ini
	ini(p,b,pgc);
	
	// Mainloop
    if(p->mpirank==0)
    cout<<"starting mainloop.SFLOW"<<endl;
    
    
//-----------MAINLOOP SFLOW----------------------------
	while(p->count<p->N45 && p->simtime<p->N41 && p->sedtime<p->S19)
	{		
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"------------------------------------"<<endl;
        cout<<p->count<<endl;
        
        cout<<"simtime: "<<setprecision(3)<<p->simtime<<endl;
		cout<<"timestep: "<<p->dt<<endl;
        
		if(p->B90>0)
		cout<<"t/T: "<<p->simtime/p->wT<<endl;
        }
        
        pflow->wavegen_2D_precalc(p,b,pgc);

        // outer loop
		double temptime=pgc->timer();
		pmom->start(p,b,pgc);
		double mtime=pgc->timer()-temptime;
    
        
        // turbulence
        pturb->start(p,b,pgc,pconvec,pdiff,psolv,pflow);
        
        // breaking waves
        int count=0;
        SLICELOOP4
        if(b->breaking(i,j)>0)
        {
        b->breaking_print(i,j)=1;
        ++count;
        }
        pgc->gcsl_start4(p,b->breaking_print,1);
    
        count=pgc->globalisum(count);
        
        if(p->mpirank==0 && (p->count%p->P12==0))
        cout<<"breaking: "<<count<<endl;
        
        if(p->A248==0)
        SLICELOOP4
        b->breaking(i,j)=0;
        
        // sediment transport
        psed->start_sflow(p,b,pgc,pflow,b->P,b->Q);
        pfsf->depth_update(p,b,pgc,b->P,b->Q,b->ws,b->eta);
        
        // 6DOF
        p6dof_sflow->start(p,pgc);

        // timesave
        pturb->ktimesave(p,b,pgc);
        pturb->etimesave(p,b,pgc);
		
        //timestep control
        p->simtime+=p->dt;
        ptime->start(p,b,pgc);
        
        
        // printer
		//print_debug(p,b,pgc);
		
		double ptime=pgc->timer();
		
        pprint->start(p,b,pgc,pflow,pturb,psed);
		pprintbed->start(p,b,pgc,psed);
		
		p->printouttime=pgc->timer()-ptime;

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
		
            if(p->mpirank==0 && (p->count%p->P12==0))
            {
            if(p->B90>0)
            cout<<"mtime: "<<setprecision(3)<<mtime<<endl;
            cout<<"wavegentime: "<<setprecision(3)<<p->wavetime<<endl;
            cout<<"printouttime: "<<setprecision(3)<<p->printouttime<<endl;
            cout<<"gctime: "<<setprecision(3)<<p->gctime<<"\t average gctime: "<<setprecision(3)<<p->gcmeantime<<endl;
            cout<<"Xtime: "<<setprecision(3)<<p->xtime<<"\t average Xtime: "<<setprecision(3)<<p->Xmeantime<<endl;		
            cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
            cout<<"timer per step: "<<setprecision(3)<<p->itertime<<endl;
            }
            
        mainlog(p);
        }
        
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavetime=0.0;
	p->lsmtime=0.0;
	p->printouttime=0.0;
    
    SLICELOOP1
    b->Pn(i,j) = b->P(i,j);
    
    SLICELOOP2
    b->Qn(i,j) = b->Q(i,j);
    
    
    if(p->umax>p->N61 || p->vmax>p->N61 || p->umax!=p->umax || p->vmax!=p->vmax)
    {
    
        if(p->mpirank==0)
        cout<<endl<<"EMERGENCY STOP  --  velocities exceeding critical value N 61"<<endl<<endl;
        
        pprint->print2D(p,b,pgc,pturb,psed);
    
    pgc->final();
    exit(0);
    }
	}

	if(p->mpirank==0)
	{
	cout<<endl<<"******************************"<<endl<<endl;

	cout<<"modelled time: "<<p->simtime<<endl;
	cout << endl;
	}
    
    pgc->final();
	
}

void sflow_f::print_debug(lexer *p, fdm2D* b, ghostcell* pgc)
{	
	
	char name[100];
	ofstream debug;
	
	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SFLOW_Log",0777);
	
	
	sprintf(name,"./REEF3D_PLS/POS-%i-%i.dat",p->count,p->mpirank+1);

	if(p->P14==0)
	sprintf(name,"/SFLOW_Debug-%i-%i.dat",p->count,p->mpirank+1);
	if(p->P14==1)
	sprintf(name,"./REEF3D_SFLOW_Log/SFLOW_Debug-%i-%i.dat",p->count,p->mpirank+1);
		
		
	debug.open(name);
	
	
	SLICELOOP4
	debug<<b->eta(i,j)<<" "<<b->bed(i,j)<<" "<<b->depth(i,j)<<endl;

	
	debug.close();
}

void sflow_f::log_ini(lexer *p)
{

	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SFLOW_Log",0777);

    if(p->mpirank==0)
    {
    if(p->P14==0)
    mainlogout.open("REEF3D_SFLOW_mainlog.dat");
    if(p->P14==1)
    mainlogout.open("./REEF3D_SFLOW_Log/REEF3D_SFLOW_mainlog.dat");

    mainlogout<<"number of cells:  "<<p->cellnumtot2D<<endl;
    mainlogout<<"#iteration \t #timestep \t #simtime \t #itertime \t #piter \t #ptime "<<endl;
    }

    
}

void sflow_f::mainlog(lexer *p)
{
	 if(p->count%p->P12==0)
	 {
     mainlogout<<fixed<<p->count<<" \t "<<setprecision(5)<<p->dt<<" \t "<<setprecision(5)<<p->simtime<<" \t ";
	 mainlogout<<fixed<<setprecision(4)<<p->itertime<<" \t ";
	 mainlogout<<p->poissoniter<<" \t "<<setprecision(4)<<p->poissontime<<" \t ";
	 mainlogout<<endl;
	 }
}

