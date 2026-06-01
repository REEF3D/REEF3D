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

#include<iomanip>
#include"nhflow_print_timeavg_wsfline.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_print_timeavg_wsfline::nhflow_print_timeavg_wsfline(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    p->Iarray(jloc,p->P146);

    maxknox=pgc->globalimax(p->knox);
    sumknox=pgc->globalisum(maxknox);

    p->Darray(xloc,p->P146+1,maxknox);
    p->Darray(wsf,p->P146+1,maxknox);
    p->Iarray(flag,p->P146+1,maxknox);
    p->Iarray(wsfpoints,p->P146+1);

    p->Darray(xloc_all,p->P146+1,sumknox);
    p->Darray(wsf_all,p->P146+1,sumknox);
    p->Iarray(flag_all,p->P146+1,sumknox);
    p->Iarray(rowflag,sumknox);

    p->Darray(yloc,p->P146);
    p->Darray(time_accum,p->P146);
    p->Darray(start_time,p->P146);
    p->Iarray(started,p->P146);
    p->Iarray(printed,p->P146);

    for(q=0;q<p->P146;++q)
    {
        time_accum[q] = 0.0;
        start_time[q] = 0.0;
        started[q] = 0;
        printed[q] = 0;
    }

    for(q=0;q<p->P146+1;++q)
    for(n=0;n<maxknox;++n)
    {
        xloc[q][n] = 1.0e20;
        wsf[q][n] = 0.0;
        flag[q][n] = 0;
    }

    for(q=0;q<p->P146+1;++q)
    for(n=0;n<sumknox;++n)
    {
        xloc_all[q][n] = 0.0;
        wsf_all[q][n] = 0.0;
        flag_all[q][n] = 0;
        rowflag[n] = 0;
    }

    for(q=0;q<p->P146;++q)
        yloc[q] = 0.0;

    location_ready = 0;

    if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_TIME-AVG-WSFLINE",0777);

}

nhflow_print_timeavg_wsfline::~nhflow_print_timeavg_wsfline()
{
    wsfout.close();
}

void nhflow_print_timeavg_wsfline::start(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, slice &f)
{
    char name[250];
    int num,check;

    if(location_ready==0)
    {
        int ready=0;

        for(q=0;q<p->P146;++q)
            if(p->simtime>=p->P146_tbegin[q] && printed[q]==0)
            {
                ready=1;
                break;
            }

        if(ready==0)
            return;

        ini_location(p,d,pgc);
        location_ready = 1;
    }

    num = p->count;

    for(q=0;q<p->P146;++q)
    {
        if(flag[q][0]<0 || printed[q]==1)
            continue;

        const double tstart = p->P146_tbegin[q];
        const double tend = p->P146_tend[q];
        const double t0 = p->simtime - p->dt_old;
        const double t1 = p->simtime;

        if(tend<=tstart)
        {
            printed[q] = 1;
            continue;
        }

        if(t1<=tstart)
            continue;

        const double effective_t0 = MAX(t0,tstart);
        const double effective_t1 = MIN(t1,tend);
        const double effective_dt = effective_t1 - effective_t0;
        if(effective_dt<=0.0)
        {
            if(t1>=tend-1.0e-12)
                printed[q] = 1;

            continue;
        }

        if(started[q]==0)
        {
            started[q] = 1;
            start_time[q] = effective_t0;
        }

        const double yline = p->P146_y[q];

        ILOOP
        {
            jloc[q] = p->posc_j(yline);

            if(p->j_dir==0)
                jloc[q] = 0;

            if(jloc[q]>=0 && jloc[q]<p->knoy)
            {
                flag[q][i] = 1;
                xloc[q][i] = p->pos_x();
                wsf[q][i] += effective_dt*(f(i,jloc[q]) + p->phimean);
            }
        }

        time_accum[q] += effective_dt;

        if(t1<tend-1.0e-12)
            continue;

        printed[q] = 1;

        if(time_accum[q]<=0.0)
            continue;

        sprintf(name,"./REEF3D_NHFLOW_TIME-AVG-WSFLINE/REEF3D-NHFLOW-TimeAvgWsfLine-%i-%i.dat",q+1,p->count);
        wsfout.open(name);

        wsfout<<"Time-Averaged Water Surface Lineprobe ID:  "<<q+1<<endl<<endl;
        wsfout<<"simtime:  "<<p->simtime<<endl;
        wsfout<<"starttime:  "<<start_time[q]<<endl;
        wsfout<<"averaging_time:  "<<time_accum[q]<<endl<<endl;
        wsfout<<"line_No     y_coord"<<endl;
        wsfout<<q+1<<"\t "<<p->Yout(0.0,yline)<<endl;

        wsfout<<endl<<endl;

        wsfout<<"X "<<q+1;
        wsfout<<"\t P "<<q+1<<endl;
        wsfout<<endl;

        wsfpoints[q]=sumknox;

            pgc->gather_double(xloc[q],maxknox,xloc_all[q],maxknox);
            pgc->gather_double(wsf[q],maxknox,wsf_all[q],maxknox);
            pgc->gather_int(flag[q],maxknox,flag_all[q],maxknox);

        if(p->mpirank==0)
        {
                sort(xloc_all[q], wsf_all[q], flag_all[q], 0, wsfpoints[q]-1);
                remove_multientry(p,xloc_all[q], wsf_all[q], flag_all[q], wsfpoints[q]);
        }

        if(p->mpirank==0)
        {
            for(n=0;n<sumknox;++n)
            {
                check=0;
                    if(flag_all[q][n]>0 && xloc_all[q][n]<1.0e20)
                {
                        wsfout<<setprecision(12)<<p->Xout(xloc_all[q][n],0.0)<<" \t ";
                        wsfout<<setprecision(12)<<(wsf_all[q][n]/time_accum[q])<<endl;
                    check=1;
                }

                if(check==1)
                    continue;
            }

            wsfout.close();
        }
    }
}

void nhflow_print_timeavg_wsfline::ini_location(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    int count;

    for(q=0;q<p->P146;++q)
    {
        count=0;
        ILOOP
        {
            if(p->j_dir==0)
                jloc[q]=0;

            if(p->j_dir==1)
                jloc[q]=p->posc_j(p->P146_y[q]);

            if(jloc[q]>=0 && jloc[q]<p->knoy)
                flag[q][count]=1;

            ++count;
        }
    }
}

void nhflow_print_timeavg_wsfline::sort(double *a, double *b, int *c, int left, int right)
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

void nhflow_print_timeavg_wsfline::remove_multientry(lexer *p, double* b, double* c, int *d, int& num)
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
