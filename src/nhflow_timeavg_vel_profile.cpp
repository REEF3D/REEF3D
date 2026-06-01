/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>
--------------------------------------------------------------------
Author: Thomas Becker
--------------------------------------------------------------------*/

#include<iomanip>
#include<cmath>
#include"nhflow_timeavg_vel_profile.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_timeavg_vel_profile::nhflow_timeavg_vel_profile(lexer *p, fdm_nhf *d) : probenum(p->P148)
{
    p->Iarray(iloc,probenum);
    p->Iarray(jloc,probenum);
    p->Iarray(flag,probenum);
    p->Darray(time_accum,probenum);
    p->Darray(start_time,probenum);
    p->Iarray(started,probenum);
    p->Iarray(printed,probenum);

    max_points = MAX(1,p->knoz);

    p->Darray(z_timeint,probenum,max_points);
    p->Darray(u_timeint,probenum,max_points);
    p->Darray(v_timeint,probenum,max_points);
    p->Darray(vel_timeint,probenum,max_points);

    for(n=0;n<probenum;++n)
    {
        time_accum[n] = 0.0;
        start_time[n] = 0.0;
        started[n] = 0;
        printed[n] = 0;

        for(q=0;q<max_points;++q)
        {
            z_timeint[n][q] = 0.0;
            u_timeint[n][q] = 0.0;
            v_timeint[n][q] = 0.0;
            vel_timeint[n][q] = 0.0;
        }
    }

    pout = new ofstream[probenum];

    if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_TIME-AVG-VEL-PROFILE",0777);

    ini_location(p,d);
}

nhflow_timeavg_vel_profile::~nhflow_timeavg_vel_profile()
{
    for(n=0;n<probenum;++n)
        pout[n].close();

    delete [] pout;
}

void nhflow_timeavg_vel_profile::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(probenum<=0)
        return;

    for(n=0;n<probenum;++n)
    {
        if(flag[n]!=1 || printed[n]==1)
            continue;

        const double tbegin = p->P148_tbegin[n];
        const double tend = p->P148_tend[n];
        const double t0 = p->simtime - p->dt_old;
        const double t1 = p->simtime;

        if(tend<=tbegin)
        {
            printed[n] = 1;
            continue;
        }

        if(t1<=tbegin)
            continue;

        const double effective_t0 = MAX(t0,tbegin);
        const double effective_t1 = MIN(t1,tend);
        const double effective_dt = effective_t1 - effective_t0;

        if(effective_dt<=0.0)
        {
            if(t1>=tend-1.0e-12)
                printed[n] = 1;

            continue;
        }

        if(started[n]==0)
        {
            started[n] = 1;
            start_time[n] = effective_t0;
        }

        i = iloc[n];
        j = jloc[n];

        for(k=0;k<p->knoz;++k)
        {
            zval = p->ZSP[IJK];
            uval = p->ccipol4V(d->U, d->WL, d->bed, p->P148_x[n], p->P148_y[n], zval);
            vval = p->ccipol4V(d->V, d->WL, d->bed, p->P148_x[n], p->P148_y[n], zval);
            velval = sqrt(uval*uval + vval*vval);

            z_timeint[n][k] += effective_dt*zval;
            u_timeint[n][k] += effective_dt*uval;
            v_timeint[n][k] += effective_dt*vval;
            vel_timeint[n][k] += effective_dt*velval;
        }

        time_accum[n] += effective_dt;

        if(t1<tend-1.0e-12)
            continue;

        printed[n] = 1;

        if(time_accum[n]<=0.0)
            continue;

        sprintf(name,"./REEF3D_NHFLOW_TIME-AVG-VEL-PROFILE/REEF3D-NHFLOW-TimeAvgVelProfile-%i-%i.dat",n+1,p->count);
        pout[n].open(name);

        pout[n]<<"TIME-AVG-VEL-Profile ID:  "<<n+1<<endl<<endl;
        pout[n]<<"x_coord     y_coord"<<endl;
        pout[n]<<n+1<<"\t "<<p->Xout(p->P148_x[n],p->P148_y[n])<<"\t "<<p->Yout(p->P148_x[n],p->P148_y[n])<<endl<<endl;
        pout[n]<<"simtime:  "<<p->simtime<<endl;
        pout[n]<<"starttime:  "<<start_time[n]<<endl;
        pout[n]<<"averaging_time:  "<<time_accum[n]<<endl<<endl;
        pout[n]<<"z\t U_avg\t V_avg\t VEL_avg"<<endl;

        for(q=0;q<p->knoz;++q)
        {
            pout[n]<<setprecision(12)
                   <<z_timeint[n][q]/time_accum[n]<<"\t "
                   <<u_timeint[n][q]/time_accum[n]<<"\t "
                   <<v_timeint[n][q]/time_accum[n]<<"\t "
                   <<vel_timeint[n][q]/time_accum[n]<<endl;
        }

        pout[n].close();
    }
}

void nhflow_timeavg_vel_profile::ini_location(lexer *p, fdm_nhf *d)
{
    int ii,jj;

    for(n=0;n<probenum;++n)
    {
        iloc[n]=p->posc_i(p->P148_x[n]);

        if(p->j_dir==0)
        {
            jloc[n]=0;
            p->P148_y[n] = 0.5*p->YP[JP];
        }

        if(p->j_dir==1)
            jloc[n]=p->posc_j(p->P148_y[n]);

        ii=iloc[n];
        jj=jloc[n];

        if(ii>=0 && ii<p->knox && jj>=0 && jj<p->knoy)
            flag[n]=1;
    }
}
