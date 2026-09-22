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
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Thomas Becker
--------------------------------------------------------------------*/

#include"nhflow_print_timeavg_wsfline_y.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<iomanip>
#include<sys/stat.h>
#include<sys/types.h>

nhflow_print_timeavg_wsfline_y::nhflow_print_timeavg_wsfline_y(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    p->Iarray(iloc,p->P147);

    maxknoy=pgc->globalimax(p->knoy);
    sumknoy=pgc->globalisum(maxknoy);

    p->Darray(yloc,p->P147,maxknoy);
    p->Darray(wsf,p->P147,maxknoy);
    p->Iarray(flag,p->P147,maxknoy);
    p->Iarray(wsfpoints,p->P147);

    p->Darray(yloc_all,p->P147,sumknoy);
    p->Darray(wsf_all,p->P147,sumknoy);
    p->Iarray(flag_all,p->P147,sumknoy);

    p->Darray(time_accum,p->P147);
    p->Darray(start_time,p->P147);
    started = new bool[p->P147];
    printed = new bool[p->P147];

    for(q=0;q<p->P147;++q)
    {
        time_accum[q] = 0.0;
        start_time[q] = 0.0;
        started[q] = false;
        printed[q] = false;
    }

    for(q=0;q<p->P147;++q)
    for(n=0;n<maxknoy;++n)
    {
        yloc[q][n] = 1.0e20;
        wsf[q][n] = 0.0;
        flag[q][n] = 0;
    }

    for(q=0;q<p->P147;++q)
    for(n=0;n<sumknoy;++n)
    {
        yloc_all[q][n] = 0.0;
        wsf_all[q][n] = 0.0;
        flag_all[q][n] = 0;
    }

    location_ready = 0;

    if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_TIME_AVG_WSFLINE_Y",0777);

}

nhflow_print_timeavg_wsfline_y::~nhflow_print_timeavg_wsfline_y()
{
    wsfout.close();
}

void nhflow_print_timeavg_wsfline_y::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, slice &f)
{
    char name[250];

    if(location_ready==0)
    {
        bool ready=false;

        for(q=0;q<p->P147;++q)
            if(p->simtime>=p->P147_tbegin[q] && printed[q]==false)
            {
                ready=true;
                break;
            }

        if(ready==false)
            return;

        ini_location(p,d,pgc);
        location_ready = 1;
    }

    for(q=0;q<p->P147;++q)
    {
        if(flag[q][0]<0 || printed[q]==true)
            continue;

        const double tbegin = p->P147_tbegin[q];
        const double tend = p->P147_tend[q];
        const double t0 = p->simtime - p->dt_old;
        const double t1 = p->simtime;

        if(tend<=tbegin)
        {
            printed[q] = true;
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
                printed[q] = true;

            continue;
        }

        if(started[q]==false)
        {
            started[q] = true;
            start_time[q] = effective_t0;
        }

        const double xline = p->P147_x[q];

        JLOOP
        {
            iloc[q] = p->posc_i(xline);

            if(iloc[q]>=0 && iloc[q]<p->knox)
            {
                flag[q][j] = 1;
                yloc[q][j] = p->pos_y();
                wsf[q][j] += effective_dt*(f(iloc[q],j) + p->phimean);
            }
        }

        time_accum[q] += effective_dt;

        if(t1<tend-1.0e-12)
            continue;

        if(time_accum[q]<=0.0)
            continue;

        wsfpoints[q]=sumknoy;

        pgc->gather_double(yloc[q],maxknoy,yloc_all[q],maxknoy);
        pgc->gather_double(wsf[q],maxknoy,wsf_all[q],maxknoy);
        pgc->gather_int(flag[q],maxknoy,flag_all[q],maxknoy);

        if(p->mpirank==0)
        {
            sort(yloc_all[q], wsf_all[q], flag_all[q], 0, wsfpoints[q]-1);
            remove_multientry(p,yloc_all[q], wsf_all[q], flag_all[q], wsfpoints[q]);

            sprintf(name,"./REEF3D_NHFLOW_TIME_AVG_WSFLINE_Y/REEF3D-NHFLOW-TimeAvgWsfLineY-%i-%i.dat",q+1,p->count);
            wsfout.open(name);

            wsfout<<"Time-Averaged Water Surface Lineprobe ID:  "<<q+1<<"\n\n"
                  <<"simtime:  "<<p->simtime<<"\n"
                  <<"starttime:  "<<start_time[q]<<"\n"
                  <<"averaging_time:  "<<time_accum[q]<<"\n\n"
                  <<"line_No     x_coord\n"
                  <<q+1<<"\t "<<p->Xout(xline,0.0)<<"\n\n"
                  <<"Y "<<q+1
                  <<"\t P "<<q+1<<"\n"<<endl;

            for(n=0;n<sumknoy;++n)
            {
                if(flag_all[q][n]>0 && yloc_all[q][n]<1.0e20)
                    wsfout<<setprecision(12)<<p->Yout(0.0,yloc_all[q][n])<<" \t "<<(wsf_all[q][n]/time_accum[q])<<endl;
            }

            wsfout.close();
            printed[q] = true;
        }
    }
}

void nhflow_print_timeavg_wsfline_y::ini_location(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    for(q=0;q<p->P147;++q)
    {
        int count=0;
        JLOOP
        {
            iloc[q]=p->posc_i(p->P147_x[q]);

            if(iloc[q]>=0 && iloc[q]<p->knox)
                flag[q][count]=1;

            ++count;
        }
    }
}

void nhflow_print_timeavg_wsfline_y::sort(double *a, double *b, int *c, int left, int right)
{
    if (left < right)
    {
        double pivot = a[right];
        int l = left;
        int r = right;

        do {
            while (a[l] < pivot) l++;
            while (a[r] > pivot) r--;

            if (l <= r) {
                double swap = a[l];
                double swapd = b[l];
                int swapc = c[l];

                a[l] = a[r];
                a[r] = swap;

                b[l] = b[r];
                b[r] = swapd;

                c[l] = c[r];
                c[r] = swapc;

                l++;
                r--;
            }
        } while (l <= r);

        sort(a,b,c, left, r);
        sort(a,b,c, l, right);
    }
}

void nhflow_print_timeavg_wsfline_y::remove_multientry(lexer *p, double* b, double* c, int *d, int& num)
{
    int oldnum=num;
    double xval=-1.12e23;
    int count=0;

    double *f,*g;
    int *h;

    p->Darray(f,num);
    p->Darray(g,num);
    p->Iarray(h,num);

    for(n=0;n<num;++n)
        g[n]=-1.12e22;

    for(n=0;n<oldnum;++n)
    {
        if(xval<=b[n]+0.001*p->DXM && xval>=b[n]-0.001*p->DXM && count>0)
            g[count-1]=MAX(g[count-1],c[n]);

        if(xval>b[n]+0.001*p->DXM || xval<b[n]-0.001*p->DXM)
        {
            f[count]=b[n];
            g[count]=c[n];
            h[count]=d[n];
            ++count;
        }

        xval=b[n];
    }

    for(n=0;n<count;++n)
    {
        b[n]=f[n];
        c[n]=g[n];
        d[n]=h[n];
    }

    p->del_Darray(f,num);
    p->del_Darray(g,num);
    p->del_Iarray(h,num);

    num=count;
}
