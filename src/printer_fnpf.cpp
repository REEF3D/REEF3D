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

#include"printer_fnpf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_force_ale.h"
#include"ioflow.h"
#include"fnpf_print_wsf.h"
#include"fnpf_print_wsf_theory.h"
#include"fnpf_print_wsfline.h"
#include"fnpf_print_wsfline_y.h"
#include"fnpf_vtp_fsf.h"
#include"fnpf_vtp_bed.h"
#include"fnpf_breaking_log.h"
#include"fnpf_print_Hs.h"
#include"fnpf_vel_probe.h"
#include"fnpf_vel_probe_theory.h"
#include"fnpf_runup.h"
#include"potentialfile_out.h"
#include"fnpf_state.h"
#include"fnpf_print_kinematics.h"
#include<sys/stat.h>
#include<sys/types.h>

printer_fnpf::printer_fnpf(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    switch (p->P10)
    {
        case 0: case 2:
            outputFormat = new vtk3D();
            break;
        case 1: default:
            outputFormat = new vtu3D();
            break;
        case 3:
            outputFormat = new vts3D();
            break;
    }

    if(p->I40==0)
    {
	p->printtime=0.0;
	p->sedprinttime=0.0;
	p->fsfprinttime=0.0;
	p->probeprinttime=0.0;
	p->stateprinttime=0.0;
    p->exportprinttime=0.0;
    }

	p->Darray(printtime_wT,p->P35);
    p->Iarray(printfsfiter_wI,p->P184);
    p->Darray(printfsftime_wT,p->P185);

	for(int qn=0; qn<p->P35; ++qn)
	printtime_wT[qn]=p->P35_ts[qn];

    for(int qn=0; qn<p->P185; ++qn)
	printfsftime_wT[qn]=p->P185_ts[qn];

    for(int qn=0; qn<p->P184; ++qn)
	printfsfiter_wI[qn]=p->P184_its[qn];


	printcount=0;

	// Create Folder
	if(p->mpirank==0)
	outputFormat->folder("FNPF");


    pwsf=new fnpf_print_wsf(p,c);

    pwsf_theory=new fnpf_print_wsf_theory(p,c,pgc);

    pwsfline=new fnpf_print_wsfline(p,c,pgc);

    pwsfline_y=new fnpf_print_wsfline_y(p,c,pgc);
    
    if(p->P65>0)
    pvel=new fnpf_vel_probe(p,c);
    
    if(p->P66>0)
    pveltheo=new fnpf_vel_probe_theory(p,c);

    if(p->P230>0)
    ppotentialfile = new potentialfile_out(p,c,pgc);

    if(p->P180==1)
	pfsf = new fnpf_vtp_fsf(p,c,pgc);

    pbed = new fnpf_vtp_bed(p,c,pgc);

    if(p->P40>0)
	pstate=new fnpf_state(p,c,pgc);

    if(p->P59==1)
    pbreaklog=new fnpf_breaking_log(p,c,pgc);
	
	if(p->P85>0)
    {
	pforce_ale = new fnpf_force_ale*[p->P85];
    
    for(n=0;n<p->P85;++n)
	pforce_ale[n]=new fnpf_force_ale(p,c,pgc,n);
    
    for(n=0;n<p->P85;++n)
    pforce_ale[n]->ini(p,c,pgc);
    }
    
    if(p->P88>0)
    {
	pkin = new fnpf_print_kinematics*[p->P88];
    
    for(n=0;n<p->P88;++n)
	pkin[n]=new fnpf_print_kinematics(p,c,pgc,n);
    
    for(n=0;n<p->P88;++n)
    pkin[n]->ini(p,c,pgc);
    }
    
    if(p->P110==1)
    phs = new fnpf_print_Hs(p,c->Hs);
    
    if(p->P140>0)
	prunup = new fnpf_runup*[p->P140];
	
	for(n=0;n<p->P140;++n)
	prunup[n]=new fnpf_runup(p,c,pgc,n);
    
    for(n=0;n<p->P140;++n)
	prunup[n]->ini(p,c,pgc);
}

printer_fnpf::~printer_fnpf()
{
}

void printer_fnpf::start(lexer* p, fdm_fnpf* c,ghostcell* pgc, ioflow *pflow)
{
    // Gages
	if(p->P51>0)
	pwsf->height_gauge(p,c,pgc,c->eta);

    if(p->P50>0)
    pwsf_theory->height_gauge(p,c,pgc,pflow);
    
    if(p->P110==1)
    phs->start(p,pgc,c->eta,c->Hs);
    
    if(p->P65>0)
	pvel->start(p,c,pgc);
    
    if(p->P66>0)
	pveltheo->start(p,c,pgc,pflow);

    // Print out based on iteration
    if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P20>0)
    {
    print_vtu(p,c,pgc);
    }

    // Print out based on time
    if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0) || (p->count==0 &&  p->P30>0.0))
    {
    print_vtu(p,c,pgc);

    p->printtime+=p->P30;
    }

    // Print out based on time interval
    if(p->P35>0)
    for(int qn=0; qn<p->P35; ++qn)
    if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
    {
    print_vtu(p,c,pgc);

    printtime_wT[qn]+=p->P35_dt[qn];
    }

    // Print FSF
    if(((p->count%p->P181==0 && p->P182<0.0 && p->P180==1 )|| (p->count==0 &&  p->P182<0.0 && p->P180==1)) && p->P181>0)
    {
    pfsf->start(p,c,pgc);
    }


    if((p->simtime>p->fsfprinttime && p->P182>0.0 && p->P180==1) || (p->count==0 &&  p->P182>0.0))
    {
    pfsf->start(p,c,pgc);
    p->fsfprinttime+=p->P182;
    }

    if(p->P180==1 && p->P184>0)
    for(int qn=0; qn<p->P184; ++qn)
    if(p->count%p->P184_dit[qn]==0 && p->count>=p->P184_its[qn] && p->count<=(p->P184_ite[qn]))
    {
    pfsf->start(p,c,pgc);
    }

    if(p->P180==1 && p->P185>0)
    for(int qn=0; qn<p->P185; ++qn)
    if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P185_ts[qn] && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
    {
    pfsf->start(p,c,pgc);

    printfsftime_wT[qn]+=p->P185_dt[qn];
    }

    // Print BED
    if(p->count==0)
    pbed->start(p,c,pgc,pflow);


    // Gages
    if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline->start(p,c,pgc,pflow,c->eta);

    if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline_y->start(p,c,pgc,pflow,c->eta);


    // Print state out based on iteration
    if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && (p->P46==0 || (p->count>=p->P46_is && p->count<<p->P46_ie)))
    {
    pstate->write(p,c,pgc);
    }

    // Print state out based on time
    if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0 && (p->P47==0 || (p->count>=p->P47_ts && p->count<<p->P47_te)))
    {
    pstate->write(p,c,pgc);

    p->stateprinttime+=p->P42;
    }

    if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
    p->probeprinttime+=p->P55;

    if(p->P59==1)
    pbreaklog->write(p,c,pgc);
	
	// ALE force    
    if(p->count>0)
    if(p->count%p->P80==0)
    for(n=0;n<p->P85;++n)
    pforce_ale[n]->start(p,c,pgc);
    
    // print kinematics    
    if(p->count>0)
    if(p->count%p->P80==0)
    for(n=0;n<p->P88;++n)
    pkin[n]->start(p,c,pgc);
    
    // Runup  
    if(p->count>0)
    for(n=0;n<p->P140;++n)
    prunup[n]->start(p,c,pgc);
}

void printer_fnpf::print_stop(lexer* p, fdm_fnpf *c, ghostcell* pgc)
{
    if(p->P180==1)
    pfsf->start(p,c,pgc);
    
    print_vtu(p,c,pgc);
}

void printer_fnpf::print_vtu(lexer* p, fdm_fnpf *c, ghostcell* pgc)
{
    if(p->P10==1||p->P10==3)
    {
        SLICELOOP4
        {
        if(c->breaking(i,j)==1)
        c->breaking_print(i,j)=1.0;

        if(c->breaking(i,j)==0)
        c->breaking_print(i,j)=0.0;
        }

        //
        pgc->start7V(p,c->Fi,c->bc,250);
        pgc->gcsl_start4(p,c->WL,50);
        pgc->gcsl_start4(p,c->bed,50);
        pgc->gcsl_start4(p,c->breaking_print,50);
        pgc->start4(p,c->test,1);

        i=-1;
        j=-1;
        if(i+p->origin_i==-1 && j+p->origin_j==-1 )
        c->WL(i,j) = c->WL(i+1,j+1);


        //----------

        outputFormat->extent(p,pgc);
        if(p->mpirank==0)
        parallel(p,pgc);

        int num=0;
        if(p->P15==1)
        num = printcount;
        if(p->P15==2)
        num = p->count;
        int rank = p->mpirank+1;
        outputFormat->fileName(name,"FNPF",num,rank);

        // Open File
        ofstream result;
        result.open(name, ios::binary);

        n=0;

        offset[n]=0;
        ++n;

        // velocity
        offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
        ++n;

        // scalars

        // Fi
        offset[n]=offset[n-1]+4*(p->pointnum)+4;
        ++n;

        // elevation
        offset[n]=offset[n-1]+4*(p->pointnum)+4;
        ++n;
        
        // test
        if(p->P23==1)
        {
        offset[n]=offset[n-1]+4*(p->pointnum)+4;
        ++n;
        }
        
        // Hs
        if(p->P110==1)
        {
        offset[n]=offset[n-1]+4*(p->pointnum)+4;
        ++n;
        }
        
        // solid
        if(p->P25==1)
        {
        offset[n]=offset[n-1]+4*(p->pointnum)+4;
        ++n;
        }

        outputFormat->offset(p,offset,n);
        //---------------------------------------------

        outputFormat->beginning(p,result);

        n=0;
        result<<"<PointData >"<<endl;
        result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        ++n;


        if(p->A10==3)
        {
        result<<"<DataArray type=\"Float32\" Name=\"Fi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        ++n;
        }

        result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        ++n;
        
        if(p->P23==1)
        {
        result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        ++n;
        }
        
        if(p->P110==1)
        {
        result<<"<DataArray type=\"Float32\" Name=\"Hs\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        ++n;
        }

        if(p->P25==1)
        {
        result<<"<DataArray type=\"Float32\" Name=\"solid\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
        ++n;
        }
        result<<"</PointData>"<<endl;


        outputFormat->ending(result,offset,n);

    //----------------------------------------------------------------------------
        result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";


    //  Velocities
        iin=3*4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
        ffn=float(c->U[FIJKp1]);

        if(k==-1 && j==-1)
        ffn=float(c->U[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof (float));


        ffn=float(c->V[FIJKp1]);

        if(k==-1 && j==-1)
        ffn=float(c->V[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof (float));


        ffn=float(c->W[FIJKp1]);

        if(k==-1 && j==-1)
        ffn=float(c->W[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof (float));
        }


    //  Fi
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        if(p->j_dir==1)
        TPLOOP
        {
        ffn=float(c->Fi[FIJKp1]);

        if(k==-1 && j==-1)
        ffn=float(c->Fi[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof (float));
        }
        
        if(p->j_dir==0)
        TPLOOP
        {
        if(j==-1)
        ffn=float(c->Fi[FIJp1Kp1]);
        
        if(j==0)
        ffn=float(c->Fi[FIJKp1]);

        if(k==-1 && j==-1)
        ffn=float(c->Fi[FIJp1Kp1]);
        result.write((char*)&ffn, sizeof (float));
        }

    //  elevation
        iin=4*(p->pointnum)*3;
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
        ffn=float(p->ZN[KP1]*c->WL(i,j) + c->bed(i,j));
        result.write((char*)&ffn, sizeof (float));
        }

    //  test
        if(p->P23==1)
        {
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
        ffn=float(p->ipol4_a(c->test));
        result.write((char*)&ffn, sizeof (float));
        }
        }
        
    //  Hs
        if(p->P110==1)
        {
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
        ffn=float(p->sl_ipol4(c->Hs));
        result.write((char*)&ffn, sizeof (float));
        }
        }

    //  solid
        if(p->P25==1)
        {
        iin=4*(p->pointnum);
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
        ffn=float(1.0);
        result.write((char*)&ffn, sizeof (float));
        }
        }

        outputFormat->structureWrite(p,c,result);

        result.close();

        ++printcount;
    }

}
