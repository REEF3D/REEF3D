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

#include"fnpf_vtu3D.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"force_ale.h"
#include"ioflow.h"
#include"fnpf_print_wsf.h"
#include"fnpf_print_wsf_theory.h"
#include"fnpf_print_wsfline.h"
#include"fnpf_print_wsfline_y.h"
#include"fnpf_vtp_fsf.h"
#include"fnpf_vtp_bed.h"
#include"fnpf_breaking_log.h"
#include"fnpf_print_Hs.h"
#include"potentialfile_out.h"
#include"fnpf_state.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_vtu3D::fnpf_vtu3D(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{	
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
	mkdir("./REEF3D_FNPF_VTU",0777);


    pwsf=new fnpf_print_wsf(p,c);

    pwsf_theory=new fnpf_print_wsf_theory(p,c,pgc);

    pwsfline=new fnpf_print_wsfline(p,c,pgc);

    pwsfline_y=new fnpf_print_wsfline_y(p,c,pgc);

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
	pforce_ale = new force_ale*[p->P85];
	
	for(n=0;n<p->P85;++n)
	pforce_ale[n]=new force_ale(p,c,pgc,n);
    
    if(p->P110==1)
    phs = new fnpf_print_Hs(p,c->Hs);

}

fnpf_vtu3D::~fnpf_vtu3D()
{
}

void fnpf_vtu3D::start(lexer* p, fdm_fnpf* c,ghostcell* pgc, ioflow *pflow)
{
    // Gages
	if(p->P51>0)
	pwsf->height_gauge(p,c,pgc,c->eta);

    if(p->P50>0)
    pwsf_theory->height_gauge(p,c,pgc,pflow);
    
    if(p->P110==1)
    phs->start(p,pgc,c->eta,c->Hs);

		// Print out based on iteration
        if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)
		{
        print_vtu(p,c,pgc);
		}

		// Print out based on time
        if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
        {
        print_vtu(p,c,pgc);

        p->printtime+=p->P30;
        }

		// Print out based on time interval
		if(p->P10==1 && p->P35>0)
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
    if((p->count==0 || p->count==p->count_statestart) && p->P85>0)
    {
    for(n=0;n<p->P85;++n)
    pforce_ale[n]->ini(p,c,pgc);
    }
    
    if(p->count>0 && p->P85>0)
    {
    for(n=0;n<p->P85;++n)
    pforce_ale[n]->start(p,c,pgc);
    }
}

void fnpf_vtu3D::print_stop(lexer* p, fdm_fnpf *c, ghostcell* pgc)
{
    if(p->P180==1)
    pfsf->start(p,c,pgc);
    
    print_vtu(p,c,pgc);
}

void fnpf_vtu3D::print_vtu(lexer* p, fdm_fnpf *c, ghostcell* pgc)
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

    pgc->dgcslpol(p,c->WL,p->dgcsl4,p->dgcsl4_count,14);
    pgc->dgcslpol(p,c->breaking_print,p->dgcsl4,p->dgcsl4_count,14);
    pgc->dgcslpol(p,c->bed,p->dgcsl4,p->dgcsl4_count,14);

    c->WL.ggcpol(p);
    c->breaking_print.ggcpol(p);

    i=-1;
    j=-1;
    if(i+p->origin_i==-1 && j+p->origin_j==-1 )
    c->WL(i,j) = c->WL(i+1,j+1);


    //----------

    if(p->mpirank==0)
    pvtu(p,pgc);

    name_iter(p,pgc);

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

	// Points
    offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
    ++n;

	// Cells
    offset[n]=offset[n-1] + 4*p->tpcellnum*8  + 4;
    ++n;
    offset[n]=offset[n-1] + 4*(p->tpcellnum)+4;
    ++n;
	offset[n]=offset[n-1] + 4*(p->tpcellnum)+4;
    ++n;
	//---------------------------------------------

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<UnstructuredGrid>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum<<"\" NumberOfCells=\""<<p->tpcellnum<<"\">"<<endl;

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


    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;

    result<<"<Cells>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Cells>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</UnstructuredGrid>"<<endl;

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

//  XYZ
	double theta_y = p->B192_1*(PI/180.0);
	double omega_y = 2.0*PI*p->B192_2;
    double waterlevel;

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    phase = omega_y*p->simtime;

	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
    waterlevel = p->sl_ipol4eta(p->wet,c->eta,c->bed)+p->wd - p->sl_ipol4(c->bed);

    zcoor = p->ZN[KP1]*waterlevel + p->sl_ipol4(c->bed);


    if(p->wet[IJ]==0)
    zcoor=c->bed(i,j);

    if(i+p->origin_i==-1 && j+p->origin_j==-1 && p->wet[(0-p->imin)*p->jmax + (0-p->jmin)]==1)
    zcoor = p->ZN[KP1]*c->WL(i,j) + c->bed(i,j);

    ffn=float( (p->XN[IP1]-p->B192_3)*cos(theta_y*sin(phase)) - (zcoor-p->B192_4)*sin(theta_y*sin(phase)) + p->B192_3);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->YN[JP1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float((p->XN[IP1]-p->B192_3)*sin(theta_y*sin(phase)) + (zcoor-p->B192_4)*cos(theta_y*sin(phase)) + p->B192_4);
	result.write((char*)&ffn, sizeof (float));
	}

//  Connectivity
    iin=4*(p->tpcellnum)*8;
    result.write((char*)&iin, sizeof (int));
    BASELOOP
    if(p->flag5[IJK]!=-20 && p->flag5[IJK]!=-30)
	{
	iin=int(c->nodeval(i-1,j-1,k-1)-1);
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i,j-1,k-1))-1;
	result.write((char*)&iin, sizeof (int));

    iin= int(c->nodeval(i,j,k-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i-1,j,k-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i-1,j-1,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i,j-1,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i,j,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(c->nodeval(i-1,j,k))-1;
	result.write((char*)&iin, sizeof (int));
	}

//  Offset of Connectivity
    iin=4*(p->tpcellnum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->tpcellnum;++n)
	{
	iin=(n+1)*8;
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*(p->tpcellnum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->tpcellnum;++n)
	{
	iin=12;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();

	++printcount;

}
