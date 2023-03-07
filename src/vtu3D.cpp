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

#include"vtu3D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"heat.h"
#include"solver.h"
#include"print_wsf.h"
#include"print_wsf_theory.h"
#include"print_wsfline_x.h"
#include"print_wsfline_y.h"
#include"force.h"
#include"vorticity_f.h"
#include"vorticity_void.h"
#include"probe_point.h"
#include"probe_pressure.h"
#include"probe_line.h"
#include"bedprobe_point.h"
#include"bedprobe_max.h"
#include"ioflow.h"
#include"data.h"
#include"concentration.h"
#include"gage_discharge_x.h"
#include"fsf_vtp.h"
#include"topo_vtp.h"
#include"state.h"
#include"bedshear_probe.h"
#include"bedshear_max.h"
#include"bedprobe_line_x.h"
#include"bedprobe_line_y.h"
#include"multiphase.h"
#include"sediment.h"
#include"sloshing_force.h"
#include"print_porous.h"
#include"export.h"
#include"flowfile_out.h"
#include<sys/stat.h>
#include<sys/types.h>

vtu3D::vtu3D(lexer* p, fdm *a, ghostcell *pgc) : nodefill(p), eta(p)
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

	pwsf=new print_wsf(p,a,pgc,0);
	pwsf_theory=new print_wsf_theory(p,a,pgc,0);
	pwsfline_x=new print_wsfline_x(p,a,pgc);
	pwsfline_y=new print_wsfline_y(p,a,pgc);
	pprobe = new probe_point(p,a,pgc);
    ppressprobe = new probe_pressure(p,a,pgc);
	pline = new probe_line(p,a,pgc);
	pq = new gage_discharge_x(p,a,pgc);

	if(p->P180==1)
	pfsf = new fsf_vtp(p,a,pgc);
    
    if(p->P190==1)
	ptopo = new topo_vtp(p,a,pgc);

    if(p->P210==1)
	pexport = new exportfile(p,a,pgc);

	if(p->P75==0)
	pvort = new vorticity_void(p,a);

	if(p->P75==1)
	pvort = new vorticity_f(p,a);

    if(p->P81>0)
	pforce = new force*[p->P81];

	if(p->P121>0)
	pbedpt = new bedprobe_point(p,a,pgc);

	if(p->P122>0)
	pbedmax = new bedprobe_max(p,a,pgc);

	if(p->P123>0)
	pbedlinex=new bedprobe_line_x(p,a,pgc);

	if(p->P124>0)
	pbedliney=new bedprobe_line_y(p,a,pgc);

	if(p->P125>0)
	pbedshear = new bedshear_probe(p,a,pgc);

	if(p->P126>0)
	pbedshearmax = new bedshear_max(p,a,pgc);

    for(n=0;n<p->P81;++n)
	pforce[n]=new force(p,a,pgc,n);

	if(p->P40>0)
	pstate=new state(p,a,pgc);

    if(p->P101>0)
	pslosh=new sloshing_force(p,a,pgc);

	if(p->B270>0 || p->B274>0 || p->B281>0 || p->B291>0 || p->B310>0 || p->B311>0)
	{
	ppor=new print_porous(p,a,pgc);
	ppor->start(p,a,pgc);
	}

    if(p->P230>0)
    pflowfile = new flowfile_out(p,a,pgc);


	p->printcount=0;

    phase=0.0;

	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_CFD_VTU",0777);
}

vtu3D::~vtu3D()
{
}

void vtu3D::ini(lexer* p, fdm* a, ghostcell* pgc)
{
}

void vtu3D::start(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
        pgc->gcparax4a(p,a->phi,5);

		// Print out based on iteration
        if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)
		{
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);
		}

		// Print out based on time
        if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
        {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

        p->printtime+=p->P30;
        }

		// Print out based on sediment time
        if((p->sedtime>p->sedprinttime && p->P34>0.0 && p->P30<0.0 && p->P10==1) || (p->count==0 &&  p->P34>0.0))
        {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

        p->sedprinttime+=p->P34;
        }

		// Print out based on time interval
		if(p->P10==1 && p->P35>0)
		for(int qn=0; qn<p->P35; ++qn)
		if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
		{
		print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

		printtime_wT[qn]+=p->P35_dt[qn];
		}


		if((p->P62>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P62>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
        pline->start(p,a,pgc,pturb);


		if(p->P50>0)
        pwsf_theory->height_gauge(p,a,pgc,pflow,a->phi);

        if(p->P51>0)
        pwsf->height_gauge(p,a,pgc,a->phi);

		if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
        pwsfline_x->wsfline(p,a,pgc,pflow);

		if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
        pwsfline_y->wsfline(p,a,pgc,pflow);

		if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
        {
        p->probeprinttime+=p->P55;
        }

		if(p->P61>0)
        pprobe->start(p,a,pgc,pturb);
        
        if(p->P64>0)
        ppressprobe->start(p,a,pgc,pturb);

		if(p->P67>0)
		pq->start(p,a,pgc);

        if((p->count==0 || p->count==p->count_statestart) && p->P81>0)
        for(n=0;n<p->P81;++n)
        pforce[n]->ini(p,a,pgc);

        if(p->count>1 && p->P81>0)
        for(n=0;n<p->P81;++n)
        pforce[n]->start(p,a,pgc);

        if(p->P101>0)
        pslosh->start(p,a,pgc);

		// sediment probes
		if(((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47) ) && p->S10>0)
		if((p->S42==1 && p->count%p->S44==0 && p->sediter%p->P120==0) || (p->S42==2 && p->simtime>=p->sedsimtime && p->sediter%p->P120==0) || (p->S42==3  && p->simtime/p->wT>=p->sedwavetime && p->sediter%p->P120==0))
		{
		if(p->P121>0)
        pbedpt->bed_gauge(p,a,pgc);

		if(p->P122>0)
        pbedmax->bed_max(p,a,pgc);

		if(p->P123>0)
        pbedlinex->start(p,a,pgc,pflow);

		if(p->P124>0)
        pbedliney->start(p,a,pgc,pflow);

		if(p->P125>0)
        pbedshear->bedshear_gauge(p,a,pgc,psed);

		if(p->P126>0)
        pbedshearmax->bedshear_maxval(p,a,pgc,psed);
		}

        // Multiphase
        pmp->print_file(p,a,pgc);

        // Print FSF
        if(((p->count%p->P181==0 && p->P182<0.0 && p->P180==1 )|| (p->count==0 &&  p->P182<0.0 && p->P180==1)) && p->P181>0)
		pfsf->start(p,a,pgc);

		if((p->simtime>p->fsfprinttime && p->P182>0.0 && p->P180==1) || (p->count==0 &&  p->P182>0.0))
        {
        pfsf->start(p,a,pgc);
        p->fsfprinttime+=p->P182;
        }

        if(p->P180==1 && p->P184>0)
		for(int qn=0; qn<p->P184; ++qn)
		if(p->count%p->P184_dit[qn]==0 && p->count>=p->P184_its[qn] && p->count<=(p->P184_ite[qn]))
		{
		pfsf->start(p,a,pgc);
		}

        if(p->P180==1 && p->P185>0)
		for(int qn=0; qn<p->P185; ++qn)
		if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P185_ts[qn] && p->simtime<=(p->P185_te[qn]+0.5*p->P185_dt[qn]))
		{
		pfsf->start(p,a,pgc);

		printfsftime_wT[qn]+=p->P185_dt[qn];
		}
        
        // Print TOPO
        if(((p->count%p->P191==0 && p->P182<0.0 && p->P190==1 )|| (p->count==0 &&  p->P192<0.0 && p->P190==1)) && p->P191>0)
		ptopo->start(p,a,pgc,psed);

		if((p->simtime>p->fsfprinttime && p->P192>0.0 && p->P190==1) || (p->count==0 &&  p->P192>0.0))
        {
        ptopo->start(p,a,pgc,psed);
        p->fsfprinttime+=p->P192;
        }

        if(p->P190==1 && p->P194>0)
		for(int qn=0; qn<p->P194; ++qn)
		if(p->count%p->P194_dit[qn]==0 && p->count>=p->P194_its[qn] && p->count<=(p->P194_ite[qn]))
		{
		ptopo->start(p,a,pgc,psed);
		}

        if(p->P190==1 && p->P195>0)
		for(int qn=0; qn<p->P195; ++qn)
		if(p->simtime>printfsftime_wT[qn] && p->simtime>=p->P195_ts[qn] && p->simtime<=(p->P195_te[qn]+0.5*p->P195_dt[qn]))
		{
		ptopo->start(p,a,pgc,psed);

		printfsftime_wT[qn]+=p->P195_dt[qn];
		}

        // Print Export
        if(p->count%p->P211==0 && p->P212<0.0 && p->P210==1)
		pexport->start(p,a,pgc);

		if((p->simtime>p->exportprinttime && p->P212>0.0 && p->P210==1) || (p->count==0 &&  p->P212>0.0))
        {
        pexport->start(p,a,pgc);
        p->exportprinttime+=p->P212;
        }

        if(p->P230>0)
        pflowfile->start(p,a,pgc,pturb);

		// Print state out based on iteration
        if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && p->P41>0)
		{
        pstate->write(p,a,pgc,pturb,psed);
		}

		// Print sate out based on time
        if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0)
        {
        pstate->write(p,a,pgc,pturb,psed);

        p->stateprinttime+=p->P42;
        }

}

void vtu3D::print_stop(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    print_vtu(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);
    
}

void vtu3D::print_vtu(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    if(p->P180==1)
	pfsf->start(p,a,pgc);
    
    print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);
}

void vtu3D::print3D(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    pgc->start4a(p,a->test,1);
    pgc->start1(p,a->u,110);
    pgc->start2(p,a->v,111);
	pgc->start3(p,a->w,112);


	pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
	pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol4(p,a->phi,14);
	pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
    //pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);

	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);
	a->press.ggcpol(p);
	a->eddyv.ggcpol(p);
	a->phi.ggcpol(p);
	a->conc.ggcpol(p);
	a->ro.ggcpol(p);
	a->visc.ggcpol(p);
	a->phi.ggcpol(p);
	a->fb.ggcpol(p);
	a->fbh4.ggcpol(p);
    //a->test.ggcpol(p);
    

    pgc->gcparacox(p,a->phi,50);
	pgc->gcparacox(p,a->phi,50);

	pgc->gcparacox(p,a->topo,150);
	pgc->gcparacox(p,a->topo,150);
    
    //pgc->start4a(p,a->topo,159);

     pgc->gcperiodicx(p,a->press,4);

    if(p->mpirank==0)
    pvtu(a,p,pgc,pturb,pheat,pdata,pconc,pmp,psed);


    name_iter(a,p,pgc);
    header(a,p,pgc);

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

		// pressure
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;

		// k and eps
	pturb->offset_vtu(p,a,pgc,result,offset,n);
		// eddyv
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
		// phi
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
		// T
	pheat->offset_vtu(p,a,pgc,result,offset,n);
    	// Multiphase
	pmp->offset_vtu(p,a,pgc,result,offset,n);
		// vorticity
	pvort->offset_vtu(p,a,pgc,result,offset,n);
		// data
	pdata->offset_vtu(p,a,pgc,result,offset,n);
		// concentration
	pconc->offset_vtu(p,a,pgc,result,offset,n);
    	// rho
	if(p->P24==1 && p->F300==0)
	{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}
        // viscosity
	if(p->P71==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}

        // omega_sig
	if(p->P72==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}

    // Fi
    if(p->A10==4)
	{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}

	if(p->P26==1)
	{
		// conc
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}
		// topo
	if(p->P27==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}
    
    	// sediment bedlaod
	if(p->P76==1)
	psed->offset_vtu_bedload(p,pgc,result,offset,n);

    	// sediment parameters 1
	if(p->P77==1)
	psed->offset_vtu_parameter1(p,pgc,result,offset,n);

    	// sediment parameters 2
	if(p->P78==1)
	psed->offset_vtu_parameter2(p,pgc,result,offset,n);

		// bed shear stress
	if(p->P79>=1)
	psed->offset_vtu_bedshear(p,pgc,result,offset,n);

    // test
    if(p->P23==1)
	{
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}
		// elevation
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;

    if(p->P25==1)
	{
		// solid
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}

	if(p->P28==1)
	{
		// floating
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}

	if(p->P29==1)
	{
		// walldist
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	}
		// end scalars

	// Points
    offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
    ++n;

	// Cells
    offset[n]=offset[n-1] + 4*p->tpcellnum*8 + 4;
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

    pturb->name_vtu(p,a,pgc,result,offset,n);

    result<<"<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;

    pheat->name_vtu(p,a,pgc,result,offset,n);
    
    pmp->name_vtu(p,a,pgc,result,offset,n);

    pvort->name_vtu(p,a,pgc,result,offset,n);

	pdata->name_vtu(p,a,pgc,result,offset,n);

	pconc->name_vtu(p,a,pgc,result,offset,n);

    if(p->P24==1 && p->F300==0)
	{
    result<<"<DataArray type=\"Float32\" Name=\"rho\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

    if(p->P71==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"viscosity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

    if(p->P72==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"omega_sig\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

    if(p->A10==4)
	{
    result<<"<DataArray type=\"Float32\" Name=\"Fi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

	if(p->P26==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"ST_conc\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

	if(p->P27==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"topo\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}
    
    if(p->P76==1)
	psed->name_vtu_bedload(p,pgc,result,offset,n);
    
    if(p->P77==1)
	psed->name_vtu_parameter1(p,pgc,result,offset,n);

    if(p->P78==1)
	psed->name_vtu_parameter2(p,pgc,result,offset,n);

	if(p->P79>=1)
	psed->name_vtu_bedshear(p,pgc,result,offset,n);

    if(p->P23==1)
	{
    result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;

    if(p->P25==1)
	{
	result<<"<DataArray type=\"Float32\" Name=\"solid\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

	if(p->P28==1)
	{
	result<<"<DataArray type=\"Float32\" Name=\"floating\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	}

	if(p->P29==1)
	{
	result<<"<DataArray type=\"Float32\" Name=\"walldist\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
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
	ffn=float(p->ipol1(a->u));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ipol2(a->v));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ipol3(a->w));
	result.write((char*)&ffn, sizeof (float));
	}

//  Pressure
	iin=4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4press(a->press));
	result.write((char*)&ffn, sizeof (float));
	}

//  turbulence
    pturb->print_3D(p,a,pgc,result);

//  eddyv
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(a->eddyv));
	result.write((char*)&ffn, sizeof (float));
	}

//  phi
	nodefill4(p,a,pgc,a->phi,eta);
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	if(p->P18==1)
	ffn=float(p->ipol4phi(a,a->phi));
	if(p->P18==2)
	ffn = float(eta(i,j,k));
	result.write((char*)&ffn, sizeof (float));
	}

//  T
    pheat->print_3D(p,a,pgc,result);
    
//  Multiphase
    pmp->print_3D(p,a,pgc,result);

//  Vorticity
    pvort->print_3D(p,a,pgc,result);

//  Data
    pdata->print_3D(p,a,pgc,result);

//  Concentration
    pconc->print_3D(p,a,pgc,result);

//  density
    if(p->P24==1 && p->F300==0)
	{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(a->ro));
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  viscosity
    if(p->P71==1)
	{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4(a->visc));
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  omega_sig
    if(p->P72==1)
	{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol3(a->omega));
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  Fi
    if(p->A10==4)
	{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4press(a->Fi));
	result.write((char*)&ffn, sizeof (float));
	}

	}

	if(p->P26==1)
	{
//  conc
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4(a->conc));
	result.write((char*)&ffn, sizeof (float));
	}
	}

	if(p->P27==1)
	{
//  topo
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
    ffn=float(p->ipol4_a(a->topo));
    //ffn = float(-a->bed(i,j)+p->ZN[KP1]);
	result.write((char*)&ffn, sizeof (float));
	}
	}
    
//  sediment bedload
	if(p->P76==1)
    psed->print_3D_bedload(p,pgc,result);
    
//  sediment parameter 1
	if(p->P77==1)
    psed->print_3D_parameter1(p,pgc,result);

//  sediment parameter 2
	if(p->P78==1)
    psed->print_3D_parameter2(p,pgc,result);

//  bed shear stress
	if(p->P79>=1)
    psed->print_3D_bedshear(p,pgc,result);

//  test
    if(p->P23==1)
	{
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(a->test));
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  elevation
	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
	ffn=float(p->pos_z()+0.5*p->DXM);
	result.write((char*)&ffn, sizeof (float));
	}

	if(p->P25==1)
	{
//  solid
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(a->solid));
	result.write((char*)&ffn, sizeof (float));
	}
	}

	if(p->P28==1 && p->X13==2)
	{
//  floating
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(a->fb));
	result.write((char*)&ffn, sizeof (float));
	}
	}

    else if (p->P28==1)
	{
//  floating
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4_a(a->fb));
	result.write((char*)&ffn, sizeof (float));
	}
	}

	if(p->P29==1)
	{
//  walldist
    iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	TPLOOP
	{
	ffn=float(p->ipol4(a->walld));
	result.write((char*)&ffn, sizeof (float));
	}
	}

//  XYZ
	double theta_y = p->B192_1*(PI/180.0);
	double omega_y = 2.0*PI*p->B192_2;

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    phase = omega_y*p->simtime;


/*
    if(p->G2==1)
    {
    pgc->gcsl_start4(p,a->WL,50);
    pgc->gcsl_start4(p,a->bed,50);

    pgc->dgcslpol(p,a->WL,p->dgcsl4,p->dgcsl4_count,14);
    pgc->dgcslpol(p,a->bed,p->dgcsl4,p->dgcsl4_count,14);

    a->WL.ggcpol(p);
    //a->test.ggcpol(p);

    i=-1;
    j=-1;
    if(i+p->origin_i==-1 && j+p->origin_j==-1 )
    a->WL(i,j) = a->WL(i+1,j+1);
    }*/


	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
        if(p->G2==0)
        zcoor=p->ZN[KP1];
/*
        if(p->G2==1)
        {
        zcoor = p->ZN[KP1]*a->WL(i,j) + a->bed(i,j);

        if(p->wet[IJ]==0 && p->flagslice4[IJ]>0)
        zcoor=a->bed(i,j);

        if(i+p->origin_i==-1 && j+p->origin_j==-1 && p->wet[(0-p->imin)*p->jmax + (0-p->jmin)]==1)
        zcoor = p->ZN[KP1]*a->WL(i,j) + a->bed(i,j);


        //cout<<"ZN: "<<p->ZN[KP1]<<" WL: "<<a->WL(i,j)<<" eta: "<<a->eta(i,j)<<" zcoor: "<<zcoor<<endl;
        }*/


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
	iin=int(a->nodeval(i-1,j-1,k-1)-1);
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval(i,j-1,k-1))-1;
	result.write((char*)&iin, sizeof (int));

    iin= int(a->nodeval(i,j,k-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval(i-1,j,k-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval(i-1,j-1,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval(i,j-1,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval(i,j,k))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(a->nodeval(i-1,j,k))-1;
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

	++p->printcount;



	pgc->start1(p,a->u,114);
    pgc->start2(p,a->v,115);
	pgc->start3(p,a->w,116);

	pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
    pgc->start4a(p,a->topo,150);

	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);

}
