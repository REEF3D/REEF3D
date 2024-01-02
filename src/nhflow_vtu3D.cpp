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

#include"nhflow_vtu3D.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"force_ale.h"
#include"ioflow.h"
#include"nhflow_print_wsf.h"
#include"nhflow_vtp_fsf.h"
#include"nhflow_vtp_bed.h"
#include"nhflow_print_wsf_theory.h"
/*#include"nhflow_print_wsfline.h"
#include"nhflow_print_wsfline_y.h"
#include"nhflow_breaking_log.h"
#include"potentialfile_out.h"
#include"nhflow_state.h"*/
#include<sys/stat.h>
#include<sys/types.h>

nhflow_vtu3D::nhflow_vtu3D(lexer* p, fdm_nhf *d, ghostcell *pgc)
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
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_NHFLOW_VTU",0777);
    
    pwsf=new nhflow_print_wsf(p,d);

    
    pwsf_theory=new nhflow_print_wsf_theory(p,d,pgc);
/*
    pwsfline=new nhflow_print_wsfline(p,d,pgc);

    pwsfline_y=new nhflow_print_wsfline_y(p,d,pgc);

    if(p->P230>0)
    ppotentialfile = new potentialfile_out(p,d,pgc);

    if(p->P40>0)
	pstate=new nhflow_state(p,d,pgc);

    if(p->P59==1)
    pbreaklog=new nhflow_breaking_log(p,d,pgc);
	
	if(p->P85>0)
	pforce_ale = new force_ale*[p->P85];
	
	for(n=0;n<p->P85;++n)
	pforce_ale[n]=new force_ale(p,d,pgc,n);*/
    
    if(p->P180==1)
	pfsf = new nhflow_vtp_fsf(p,d,pgc);

    pbed = new nhflow_vtp_bed(p,d,pgc);

}

nhflow_vtu3D::~nhflow_vtu3D()
{
}

void nhflow_vtu3D::start(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow *pflow)
{
    // Gages
	if(p->P51>0)
	pwsf->height_gauge(p,d,pgc,d->eta);

    if(p->P50>0)
    pwsf_theory->height_gauge(p,d,pgc,pflow);

		// Print out based on iteration
        if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)
		{
        print_vtu(p,d,pgc);
		}

		// Print out based on time
        if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
        {
        print_vtu(p,d,pgc);

        p->printtime+=p->P30;
        }

		// Print out based on time interval
		if(p->P10==1 && p->P35>0)
		for(int qn=0; qn<p->P35; ++qn)
		if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
		{
		print_vtu(p,d,pgc);

		printtime_wT[qn]+=p->P35_dt[qn];
		}

        // Print FSF
		if(((p->count%p->P181==0 && p->P182<0.0 && p->P180==1 )|| (p->count==0 &&  p->P182<0.0 && p->P180==1)) && p->P181>0)
        {
		pfsf->start(p,d,pgc);
        }


		if((p->simtime>p->fsfprinttime && p->P182>0.0 && p->P180==1) || (p->count==0 &&  p->P182>0.0))
        {
        pfsf->start(p,d,pgc);
        p->fsfprinttime+=p->P182;
        }

        if(p->P180==1 && p->P184>0)
		for(int qn=0; qn<p->P184; ++qn)
		if(p->count%p->P184_dit[qn]==0 && p->count>=p->P184_its[qn] && p->count<=(p->P184_ite[qn]))
		{
		pfsf->start(p,d,pgc);
		}

         if(p->P180==1 && p->P185>0)
		for(int qn=0; qn<p->P185; ++qn)
		if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P185_ts[qn] && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
		{
		pfsf->start(p,d,pgc);

		printfsftime_wT[qn]+=p->P185_dt[qn];
		}

        // Print BED
        if(p->count==0)
		pbed->start(p,d,pgc);

/*
    // Gages
    if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline->start(p,d,pgc,pflow,d->eta);

    if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline_y->start(p,d,pgc,pflow,d->eta);


    // Print state out based on iteration
    if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 &&)
    {
    pstate->write(p,d,pgc);
    }

    // Print sate out based on time
    if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0)
    {
    pstate->write(p,d,pgc);

    p->stateprinttime+=p->P42;
    }

    if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
    p->probeprinttime+=p->P55;

    if(p->P59==1)
    pbreaklog->write(p,d,pgc);
	
	// ALE force
	  if((p->count==0 || p->count==p->count_statestart) && p->P85>0)
	  {
		for(n=0;n<p->P85;++n)
        pforce_ale[n]->ini(p,d,pgc);
	  }
        if(p->count>0 && p->P85>0)
		{
        for(n=0;n<p->P85;++n)
        pforce_ale[n]->start(p,d,pgc);
		}
        */
}

void nhflow_vtu3D::print_stop(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow *pflow)
{
    print_vtu(p,d,pgc);
    
    if(p->P180==1)
    pfsf->start(p,d,pgc);
    
}

void nhflow_vtu3D::print_vtu(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
    /*
    - U, V, W
    - P
    - test
    - breaking
    */
    
    SLICELOOP4
    {
    if(d->breaking(i,j)==1)
    d->breaking_print(i,j)=1.0;

    if(d->breaking(i,j)==0)
    d->breaking_print(i,j)=0.0;
    }
    
    //
    //pgc->gcsl_start4(p,d->WL,50);
    pgc->gcsl_start4(p,d->bed,50);
    pgc->gcsl_start4(p,d->breaking_print,50);
    pgc->start4V(p,d->test,50);
    //pgc->start4(p,d->test,1);


    pgc->dgcslpol(p,d->WL,p->dgcsl4,p->dgcsl4_count,14);
    pgc->dgcslpol(p,d->breaking_print,p->dgcsl4,p->dgcsl4_count,14);
    pgc->dgcslpol(p,d->bed,p->dgcsl4,p->dgcsl4_count,14);

    d->WL.ggcpol(p);
    d->breaking_print.ggcpol(p);

    i=-1;
    j=-1;
    if(i+p->origin_i==-1 && j+p->origin_j==-1 )
    d->WL(i,j) = d->WL(i+1,j+1);


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

    // P
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
    
    // omega_sig
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


    result<<"<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
    result<<"<DataArray type=\"Float32\" Name=\"omega_sig\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;


    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    
    if(p->P23==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
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
    if(p->j_dir==0)
    {
    jj=j;
    j=0;
	ffn=float(0.5*(d->U[IJK]+d->U[IJKp1]));
    j=jj;
    }
    
    
    if(p->j_dir==1)
	ffn=float(0.25*(d->U[IJK]+d->U[IJKp1]+d->U[IJp1K]+d->U[IJp1Kp1]));
    
	result.write((char*)&ffn, sizeof (float));


	if(p->j_dir==0)
    {
    jj=j;
    j=0;
	ffn=float(0.5*(d->V[IJK]+d->V[IJKp1]));
    j=jj;
    }

    if(p->j_dir==1)
	ffn=float(0.25*(d->V[IJK]+d->V[IJKp1]+d->V[IJp1K]+d->V[IJp1Kp1]));
    
	result.write((char*)&ffn, sizeof (float));


	if(p->j_dir==0)
    {
    jj=j;
    j=0;
	ffn=float(0.5*(d->W[IJK]+d->W[Im1JK]));
    j=jj;
    }

    if(p->j_dir==1)
	ffn=float(0.25*(d->W[IJK]+d->W[Im1JK]+d->W[IJK]+d->W[IJm1K]));
    
	result.write((char*)&ffn, sizeof (float));
	}

//  P
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
        
    if(p->A520<10)
    {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(d->P[FIJKp1]);
        j=jj;
        }

        if(p->j_dir==1)
        ffn=float(0.5*(d->P[FIJKp1]+d->P[FIJKp1]));
    }
    
    if(p->A520>=10)
    {
        if(p->j_dir==0)
        {
        jj=j;
        j=0;
        ffn=float(0.25*(d->P[IJK]+d->P[IJKp1]+d->P[Ip1JK]+d->P[Ip1JKp1]));
        j=jj;
        }

        if(p->j_dir==1)
        ffn=float(0.25*(d->P[IJK]+d->P[IJKp1]+d->P[IJp1K]+d->P[IJp1Kp1]));
    }
    
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Omega_sig
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
    if(p->j_dir==0)
    {
    jj=j;
    j=0;
	ffn=float((d->omegaF[FIJKp1]));
    j=jj;
    }

    if(p->j_dir==1)
	ffn=float(0.5*(d->omegaF[FIJKp1]+d->omegaF[FIJp1Kp1]));
    
	result.write((char*)&ffn, sizeof (float));
	}

//  elevation
	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
	ffn=float(p->ZN[KP1]*d->WL(i,j) + d->bed(i,j));
	result.write((char*)&ffn, sizeof (float));
	}

//  test
    if(p->P23==1)
	{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	if(p->j_dir==0)
    {
    jj=j;
    j=0;
	ffn=float(0.5*(d->test[IJK]+d->test[IJKp1]));
    j=jj;
    }
    
    if(p->j_dir==1)
	ffn=float(0.25*(d->test[IJK]+d->test[IJKp1]+d->test[IJp1K]+d->test[IJp1Kp1]));
    
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  XYZ
	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
    zcoor = p->ZN[KP1]*p->sl_ipol4(d->WL) + p->sl_ipol4(d->bed);

    if(p->wet[IJ]==0)
    zcoor=p->sl_ipol4(d->bed);
    
    if(i+p->origin_i==-1 && j+p->origin_j==-1 && p->wet[(0-p->imin)*p->jmax + (0-p->jmin)]==1)
    zcoor = p->ZN[KP1]*d->WL(i,j) + d->bed(i,j);

    // -- 
    ffn=float(p->XN[IP1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->YN[JP1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(zcoor);
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Connectivity
    iin=4*(p->tpcellnum)*8;
    result.write((char*)&iin, sizeof (int));
    BASELOOP
    if(p->flag5[IJK]!=-20 && p->flag5[IJK]!=-30)
	{
    iin=int(d->NODEVAL[Im1Jm1Km1])-1;
	result.write((char*)&iin, sizeof (int));

    iin=int(d->NODEVAL[IJm1Km1])-1;
	result.write((char*)&iin, sizeof (int));

    iin= int(d->NODEVAL[IJKm1])-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->NODEVAL[Im1JKm1])-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->NODEVAL[Im1Jm1K])-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->NODEVAL[IJm1K])-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->NODEVAL[IJK])-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(d->NODEVAL[Im1JK])-1;
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
