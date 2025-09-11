/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"printer_CFD.h"
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
#include"ioflow.h"
#include"data.h"
#include"concentration.h"
#include"gage_discharge_x.h"
#include"gage_discharge_window_x.h"
#include"fsf_vtp.h"
#include"topo_vtp.h"
#include"cfd_state.h"

#include"probe_vel.h"
#include"probe_vel_theory.h"
#include"multiphase.h"
#include"sediment.h"
#include"sloshing_force.h"
#include"print_porous.h"
#include"flowfile_out.h"
#include"print_averaging_f.h"
#include"print_averaging_v.h"
#include<sys/stat.h>
#include<sys/types.h>

printer_CFD::printer_CFD(lexer* p, fdm *a, ghostcell *pgc)
{
    switch (p->P10)
    {
        case vtk3D::type::none:
            outputFormat = new vtk3D();
            break;
        case vtk3D::type::vtu: default:
            outputFormat = new vtu3D();
            break;
        case vtk3D::type::vtr:
            outputFormat = new vtr3D();
            break;
        case vtk3D::type::vts:
            outputFormat = new vts3D();
            break;
    }

    if(p->I40==0)
    {
        p->printtime=0.0;
        p->sedprinttime=0.0;
        p->fsfprinttime=0.0;
        p->fsfsedprinttime=0.0;
        p->probeprinttime=0.0;
        p->stateprinttime=0.0;
        p->exportprinttime=0.0;
    }

    p->printcount=0;

    p->Darray(printtime_wT,p->P35);
    p->Darray(printfsftime_wT,p->P185);

    for(int qn=0; qn<p->P35; ++qn)
        printtime_wT[qn] = p->P35_ts[qn];

    for(int qn=0; qn<p->P185; ++qn)
        printfsftime_wT[qn] = p->P185_ts[qn];

    pwsf = new print_wsf(p,a,pgc,0);
    pwsf_theory = new print_wsf_theory(p,a,pgc,0);
    pwsfline_x = new print_wsfline_x(p,a,pgc);
    pwsfline_y = new print_wsfline_y(p,a,pgc);
    pprobe = new probe_point(p,a,pgc);
    ppressprobe = new probe_pressure(p,a,pgc);
    pline = new probe_line(p,a,pgc);
    pq = new gage_discharge_x(p,a,pgc);
    pqw = new gage_discharge_window_x(p,a,pgc);

    if(p->P21==0)
        pmean = new print_averaging_v(p,a,pgc);
    else if(p->P21==1)
        pmean = new print_averaging_f(p,a,pgc);

    if(p->P180==1)
        pfsf = new fsf_vtp(p,a,pgc);

    if(p->P190==1)
        ptopo = new topo_vtp(p,a,pgc);

    if(p->P65>0)
        pvel = new probe_vel(p,a);

    if(p->P66>0)
        pveltheo = new probe_vel_theory(p,a);

    if(p->P75==0)
        pvort = new vorticity_void(p,a);

    if(p->P75==1)
        pvort = new vorticity_f(p,a);

    if(p->P81>0)
    {
        pforce = new force*[p->P81];

        for(n=0;n<p->P81;++n)
            pforce[n] = new force(p,a,pgc,n);
    }

    // if(p->P37>0)
    //     pstate_restart = new cfd_state(p,a,pgc,0);

    if(p->P40>0)
        pstate = new cfd_state(p,a,pgc,0);

    if(p->P101>0)
        pslosh = new sloshing_force(p,a,pgc);

    if(p->B270>0 || p->B274>0 || p->B281>0 || p->B282>0 || p->B291>0 || p->B310>0 || p->B321>0 || p->B322>0 || p->B311>0)
    {
        ppor = new print_porous(p,a,pgc);
        ppor->start(p,a,pgc);
    }

    if(p->P230>0)
        pflowfile = new flowfile_out(p,a,pgc);

    // Create Folder
    if(p->mpirank==0)
        outputFormat->folder("CFD");
}

void printer_CFD::start(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    pgc->gcparax4a(p,a->phi,5);

    pmean->averaging(p,a,pgc,pheat);

    // Print out based on iteration
    if(p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P20>0)
    {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);
    }

    // Print out based on time
    if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0) || (p->count==0 &&  p->P30>0.0))
    {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

        p->printtime+=p->P30;
    }

    // Print out based on sediment time
    if((p->sedtime>p->sedprinttime && p->P34>0.0 && p->P30<0.0) || (p->count==0 &&  p->P34>0.0))
    {
        print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);

        p->sedprinttime+=p->P34;
    }

    // Print out based on time interval
    if(p->P35>0)
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

    if(p->P51>0 && p->F80!=4)
        pwsf->height_gauge(p,a,pgc,a->phi);

    if(p->P51>0 && p->F80==4)
        pwsf->height_gauge(p,a,pgc,a->vof);

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

    if(p->P65>0)
        pvel->start(p,a,pgc);

    if(p->P66>0)
        pveltheo->start(p,a,pgc,pflow);

    if(p->P167>0)
        pq->start(p,a,pgc);

    if(p->P168>0)
        pqw->start(p,a,pgc);

    if((p->count==0 || p->count==p->count_statestart) && p->P81>0)
        for(n=0;n<p->P81;++n)
            pforce[n]->ini(p,a,pgc);

    if(p->count>1 && p->P81>0)
        for(n=0;n<p->P81;++n)
            pforce[n]->start(p,a,pgc);

    if(p->P101>0)
        pslosh->start(p,a,pgc);

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

    if(p->P230>0)
        pflowfile->start(p,a,pgc,pturb);

    // Print state out based on iteration
    if(p->count%p->P41==0 && p->P42<0.0 && p->P40>0 && p->P41>0 && (p->P46==0 || (p->count>=p->P46_is && p->count<<p->P46_ie)))
    {
        pstate->write(p,a,pgc,pturb,psed);
    }

    // Print state out based on time
    if((p->simtime>p->stateprinttime && p->P42>0.0 || (p->count==0 &&  p->P42>0.0)) && p->P40>0 && (p->P47==0 || (p->count>=p->P47_ts && p->count<<p->P47_te)))
    {
        pstate->write(p,a,pgc,pturb,psed);

        p->stateprinttime+=p->P42;
    }

    /* // Print state restart out based on iteration
    if(p->count%p->P38==0 && p->P37>0)
    {
    pstate_restart->write(p,a,pgc,pturb,psed);
    }*/

}

void printer_CFD::print_stop(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, ioflow *pflow, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    if(p->P180==1)
        pfsf->start(p,a,pgc);

    print3D(a,p,pgc,pturb,pheat,psolv,pdata,pconc,pmp,psed);
}

void printer_CFD::print3D(fdm* a,lexer* p,ghostcell* pgc, turbulence *pturb, heat *pheat, solver *psolv, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    if(p->P10!=vtk3D::type::none)
    {
        pgc->start4a(p,a->test,1);
        pgc->start1(p,a->u,110);
        pgc->start2(p,a->v,111);
        pgc->start3(p,a->w,112);

        pgc->gcparacox(p,a->phi,50);
        pgc->gcparacox(p,a->phi,50);

        pgc->gcparacox(p,a->topo,150);
        pgc->gcparacox(p,a->topo,150);

        //pgc->start4a(p,a->topo,159);

        pgc->gcperiodicx(p,a->press,4);

        outputFormat->extent(p,pgc);
        if(p->mpirank==0)
                parallel(a,p,pgc,pturb,pheat,pdata,pconc,pmp,psed);

        int num=0;
        if(p->P15==1)
            num = p->printcount;
        else if(p->P15==2)
            num = p->count;
        outputFormat->fileName(name,sizeof(name),"CFD",num,p->mpirank+1);

        // Open File
        ofstream result;
        result.open(name, ios::binary);

        n=0;

        offset[n]=0;
        ++n;

        // velocity
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum*3+sizeof(int);
        ++n;

        pmean->offset_ParaView(p,offset,n);

        // scalars

        // pressure
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
        ++n;
        // k and eps
        pturb->offset_ParaView(p,offset,n);
        // eddyv
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
        ++n;
        // phi
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
        ++n;
        // T
        pheat->offset_ParaView(p,offset,n);
        // Multiphase
        pmp->offset_ParaView(p,offset,n);
        // vorticity
        pvort->offset_ParaView(p,offset,n);
        // data
        pdata->offset_ParaView(p,offset,n);
        // concentration
        pconc->offset_ParaView(p,offset,n);
        // rho
        if(p->P24==1 && p->F300==0)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // viscosity
        if(p->P71==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // VOF
        if(p->P72==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // Fi
        if(p->A10==4)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // conc
        if(p->P26==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // topo
        if(p->P27==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // sediment bedlaod
        if(p->P76==1)
            psed->offset_ParaView_bedload(p,offset,n);
        // sediment parameters 1
        if(p->P77==1)
            psed->offset_ParaView_parameter1(p,offset,n);
        // sediment parameters 2
        if(p->P78==1)
            psed->offset_ParaView_parameter2(p,offset,n);
        // bed shear stress
        if(p->P79>=1)
            psed->offset_ParaView_bedshear(p,offset,n);
        // test
        if(p->P23==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // elevation
        offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
        ++n;
        // solid
        if(p->P25==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // floating
        if(p->P28==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
            ++n;
        }
        // walldist
        if(p->P29==1)
        {
            offset[n]=offset[n-1]+sizeof(float)*p->pointnum+sizeof(int);
        ++n;
        }
        // end point data
        // VOF_C
        if(p->P72==1)
        {
            offset[n]=offset[n-1]+4*(p->cellnum)+4;
            ++n;
        }
        // end scalars

        // Format specific structure
        outputFormat->offset(p,offset,n);
        //---------------------------------------------

        outputFormat->beginning(p,result);

        n=0;
        result<<"<PointData>\n";
        result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;

        pmean->name_ParaView(p,a,pgc,result,offset,n);

        result<<"<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;

        pturb->name_ParaView(p,a,pgc,result,offset,n);

        result<<"<DataArray type=\"Float32\" Name=\"eddyv\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
        result<<"<DataArray type=\"Float32\" Name=\"phi\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;

        pheat->name_ParaView(p,a,pgc,result,offset,n);

        pmp->name_ParaView(p,a,pgc,result,offset,n);

        pvort->name_ParaView(p,a,pgc,result,offset,n);

        pdata->name_ParaView(p,a,pgc,result,offset,n);

        pconc->name_ParaView(p,a,pgc,result,offset,n);

        if(p->P24==1 && p->F300==0)
        {
            result<<"<DataArray type=\"Float32\" Name=\"rho\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P71==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"viscosity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P72==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"VOF\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->A10==4)
        {
            result<<"<DataArray type=\"Float32\" Name=\"Fi\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P26==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"ST_conc\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P27==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"topo\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P76==1)
            psed->name_ParaView_bedload(p,pgc,result,offset,n);

        if(p->P77==1)
            psed->name_ParaView_parameter1(p,pgc,result,offset,n);

        if(p->P78==1)
            psed->name_ParaView_parameter2(p,pgc,result,offset,n);

        if(p->P79>=1)
            psed->name_ParaView_bedshear(p,pgc,result,offset,n);

        if(p->P23==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        result<<"<DataArray type=\"Float32\" Name=\"elevation\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;

        if(p->P25==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"solid\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P28==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"floating\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }

        if(p->P29==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"walldist\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
            ++n;
        }
        result<<"</PointData>\n";
        result<<"<CellData>\n";
        if(p->P72==1)
        {
            result<<"<DataArray type=\"Float32\" Name=\"VOF_C\"  format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
            ++n;
        }
        result<<"</CellData>\n";

        outputFormat->ending(result,offset,n);

        //----------------------------------------------------------------------------


        //  Velocities
        iin=3*sizeof(int)*p->pointnum;
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

        //  time average flow parameters
        pmean->print_3D(p,a,pgc,result);

        //  Pressure
        iin=sizeof(float)*p->pointnum;
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->ipol4press(a->press)-p->pressgage);
            result.write((char*)&ffn, sizeof (float));
        }

        //  turbulence
        pturb->print_3D(p,a,pgc,result);

        //  eddyv
        iin=sizeof(float)*p->pointnum;
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->ipol4_a(a->eddyv));
            result.write((char*)&ffn, sizeof (float));
        }

        //  phi
        iin=sizeof(float)*p->pointnum;
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->ipol4phi(a,a->phi));
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
            iin=sizeof(float)*p->pointnum;
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
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4(a->visc));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  VOF
        if(p->P72==1)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4(a->vof));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  Fi
        if(p->A10==4)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4press(a->Fi));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  conc
        if(p->P26==1)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4(a->conc));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  topo
        if(p->P27==1)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->topo));
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
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->test));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  elevation
        iin=sizeof(float)*p->pointnum*3;
        result.write((char*)&iin, sizeof (int));
        TPLOOP
        {
            ffn=float(p->pos_z()+0.5*p->DZN[KP]);
            result.write((char*)&ffn, sizeof (float));
        }

        //  solid
        if(p->P25==1)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->solid));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  floating
        if(p->P28==1)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->fb));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  walldist
        if(p->P29==1)
        {
            iin=sizeof(float)*p->pointnum;
            result.write((char*)&iin, sizeof (int));
            TPLOOP
            {
                ffn=float(p->ipol4_a(a->walld));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        //  VOF_C
        if(p->P72==1)
        {
            iin=4*(p->cellnum);
            result.write((char*)&iin, sizeof (int));
            BASEREVLOOP
            {
                ffn=float(a->vof(i,j,k));
                result.write((char*)&ffn, sizeof (float));
            }
        }

        outputFormat->structureWrite(p,a,result);

        result.close();

        ++p->printcount;

        pgc->start1(p,a->u,114);
        pgc->start2(p,a->v,115);
        pgc->start3(p,a->w,116);

        pgc->start4a(p,a->topo,150);
    }
}
