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

            cout<<"simtime: "<<p->simtime<<endl;
            cout<<"timestep: "<<p->dt<<endl;

            if(p->B90>0 && p->B92<=11)
            cout<<"t/T: "<<p->simtime/p->wT<<endl;
            
            if(p->B90>0 && p->B92>11)
            cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }

        pflow->flowfile(p,a,pgc,pturb);

        pflow->wavegen_precalc(p,pgc);

        fill_vel(p,a,pgc);
        
        // Benchmark cases
        pbench->start(p,a,pgc,pconvec);
/*        
        // Free-surface computation
        field1 uold(p); ULOOP uold(i,j,k) = a->u(i,j,k); pgc->start1(p,uold,10);
        field2 vold(p); VLOOP vold(i,j,k) = a->v(i,j,k); pgc->start2(p,vold,11);
        field3 wold(p); WLOOP wold(i,j,k) = a->w(i,j,k); pgc->start3(p,wold,12);

        double nx, ny, nz, norm, dist, xloc, yloc, zloc;
        double uvel_s,vvel_s,wvel_s,un_s,un_x_s,un_y_s,un_z_s,uvel_f,vvel_f,wvel_f,un_f,ut_x_f,ut_y_f,ut_z_f;

        ULOOP
        {
            if (a->fbh1(i,j,k) > 0.0)
            {
                nx = -(a->fbh1(i+1,j,k) - a->fbh1(i-1,j,k))/(2.0*p->DXP[IP]);
                ny = -(a->fbh1(i,j+1,k) - a->fbh1(i,j-1,k))/(2.0*p->DYN[JP]);
                nz = -(a->fbh1(i,j,k+1) - a->fbh1(i,j,k-1))/(2.0*p->DZN[KP]);

                norm = sqrt(nx*nx + ny*ny + nz*nz);
                
                nx /= norm > 1.0e-20 ? norm : 1.0e20;
                ny /= norm > 1.0e-20 ? norm : 1.0e20;
                nz /= norm > 1.0e-20 ? norm : 1.0e20;
                
                dist = 1.0*p->X41*1.0/3.0*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
                
                xloc = p->pos1_x() + nx*dist;
                yloc = p->pos1_y() + ny*dist;
                zloc = p->pos1_z() + nz*dist;
  */              
                /*
                uvel_s = p->ccipol1_a(uold,p->pos1_x(),p->pos1_y(),p->pos1_z());
                vvel_s = p->ccipol2_a(vold,p->pos2_x(),p->pos2_y(),p->pos2_z());
                wvel_s = p->ccipol3_a(wold,p->pos3_x(),p->pos3_y(),p->pos3_z());

                un_s = uvel_s*nx + vvel_s*ny + wvel_s*nz; 
                un_x_s = un_s*nx;
                un_y_s = un_s*ny;
                un_z_s = un_s*nz;
                
                xloc = p->pos1_x() + nx*dist;
                yloc = p->pos1_y() + ny*dist;
                zloc = p->pos1_z() + nz*dist;
            
                uvel_f = p->ccipol1_a(uold,xloc,yloc,zloc);
                vvel_f = p->ccipol2_a(vold,xloc,yloc,zloc);
                wvel_f = p->ccipol3_a(wold,xloc,yloc,zloc);

                un_f = uvel_f*nx + vvel_f*ny + wvel_f*nz; 
                ut_x_f = uvel_f - un_f*nx;
                ut_y_f = vvel_f - un_f*ny;
                ut_z_f = uvel_f - un_f*nz;

                a->u(i,j,k) = un_x_s + ut_x_f; 
                */ 
    /*            
                a->u(i,j,k) = p->ccipol1_a(uold,xloc,yloc,zloc);
            }
        }

        VLOOP
        {
            if (a->fbh2(i,j,k) > 0.0)
            {
                nx = -(a->fbh2(i+1,j,k) - a->fbh2(i-1,j,k))/(2.0*p->DXN[IP]);
                ny = -(a->fbh2(i,j+1,k) - a->fbh2(i,j-1,k))/(2.0*p->DYP[JP]);
                nz = -(a->fbh2(i,j,k+1) - a->fbh2(i,j,k-1))/(2.0*p->DZN[KP]);

                norm = sqrt(nx*nx + ny*ny + nz*nz);
                
                nx /= norm > 1.0e-20 ? norm : 1.0e20;
                ny /= norm > 1.0e-20 ? norm : 1.0e20;
                nz /= norm > 1.0e-20 ? norm : 1.0e20;
                
                dist = p->X41*1.0/3.0*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
                
                xloc = p->pos2_x() + nx*dist;
                yloc = p->pos2_y() + ny*dist;
                zloc = p->pos2_z() + nz*dist;
                
                a->v(i,j,k) = p->ccipol2_a(vold,xloc,yloc,zloc);
            }
        }
        
        WLOOP
        {
            if (a->fbh3(i,j,k) > 0.0)
            {
                nx = -(a->fbh3(i+1,j,k) - a->fbh3(i-1,j,k))/(2.0*p->DXN[IP]);
                ny = -(a->fbh3(i,j+1,k) - a->fbh3(i,j-1,k))/(2.0*p->DYN[JP]);
                nz = -(a->fbh3(i,j,k+1) - a->fbh3(i,j,k-1))/(2.0*p->DZP[KP]);

                norm = sqrt(nx*nx + ny*ny + nz*nz);
                
                nx /= norm > 1.0e-20 ? norm : 1.0e20;
                ny /= norm > 1.0e-20 ? norm : 1.0e20;
                nz /= norm > 1.0e-20 ? norm : 1.0e20;
                
                dist = 1.0*p->X41*1.0/3.0*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);

                xloc = p->pos3_x() + nx*dist;
                yloc = p->pos3_y() + ny*dist;
                zloc = p->pos3_z() + nz*dist;
                
                a->w(i,j,k) = p->ccipol3_a(wold,xloc,yloc,zloc);
            }
        }

        pgc->start1(p,a->u,10);
        pgc->start2(p,a->v,11);
        pgc->start3(p,a->w,12); 

        pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,ppart,a->phi);
        poneph->update(p,a,pgc,pflow);
       
        ULOOP a->u(i,j,k) = uold(i,j,k); pgc->start1(p,a->u,10);
        VLOOP a->v(i,j,k) = vold(i,j,k); pgc->start2(p,a->v,11);
        WLOOP a->w(i,j,k) = wold(i,j,k); pgc->start3(p,a->w,12);
*/

        pfsf->start(a,p, pfsfdisc,psolv,pgc,pflow,preini,ppart,a->phi);
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
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;

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

