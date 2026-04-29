/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"nhflow_depavg_vel_lineprobe.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>
#include<iomanip>
#include<sys/stat.h>
#include<sys/types.h>

nhflow_depavg_vel_lineprobe::nhflow_depavg_vel_lineprobe(lexer *p, fdm_nhf *d) : probenum(p->P144)
{
    p->Iarray(flag,probenum);
    p->Darray(line_length,probenum);
    p->Darray(normal_x,probenum);
    p->Darray(normal_y,probenum);
    p->Iarray(point_count,probenum);

    max_points = 1;
    for(n=0;n<probenum;++n)
    {
        point_count[n] = MAX(1,p->P144_n[n]);
        max_points = MAX(max_points,point_count[n]);
    }

    p->Darray(xpt,probenum,max_points);
    p->Darray(ypt,probenum,max_points);
    p->Darray(dist,probenum,max_points);
    p->Iarray(iloc,probenum,max_points);
    p->Iarray(jloc,probenum,max_points);

    if(p->mpirank==0)
        mkdir("./REEF3D_NHFLOW_DEP-AVG-VEL-LINE",0777);

    pout = new ofstream[probenum];

    ini_location(p,d);
}

nhflow_depavg_vel_lineprobe::~nhflow_depavg_vel_lineprobe()
{
    for(n=0;n<probenum;++n)
        pout[n].close();

    delete [] pout;
}

void nhflow_depavg_vel_lineprobe::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(probenum<=0)
        return;

    for(n=0;n<probenum;++n)
    {
        if(flag[n]!=1)
            continue;

        sprintf(name,"./REEF3D_NHFLOW_DEP-AVG-VEL-LINE/REEF3D-NHFLOW-DepAvgVel-LineProbe-%i-%i.dat",n+1,p->count);
        pout[n].open(name);

        pout[n]<<"Depth-Averaged Velocity Lineprobe ID:  "<<n+1<<endl<<endl;
        pout[n]<<"simtime:  "<<p->simtime<<endl<<endl;
        pout[n]<<"dist_along_the_line x_point y_point u_depth_avg v_depth_avg norm_vel_depth_avg"<<endl;

        for(q=0;q<point_count[n];++q)
        {
            double uavg_local = -1.0e20;
            double vavg_local = -1.0e20;

            if(iloc[n][q]>=0 && jloc[n][q]>=0)
            {
                i = iloc[n][q];
                j = jloc[n][q];

                double usum = 0.0;
                double vsum = 0.0;
                double wsum = 0.0;

                for(k=0;k<p->knoz;++k)
                {
                    const double weight = p->DZN[KP];
                    const double zpos = p->ZSP[IJK];

                    const double uval = p->ccipol4V(d->U, d->WL, d->bed, xpt[n][q], ypt[n][q], zpos);
                    const double vval = p->ccipol4V(d->V, d->WL, d->bed, xpt[n][q], ypt[n][q], zpos);

                    usum += weight*uval;
                    vsum += weight*vval;
                    wsum += weight;
                }

                if(wsum>1.0e-12)
                {
                    uavg_local = usum/wsum;
                    vavg_local = vsum/wsum;
                }
            }

            const double uavg = pgc->globalmax(uavg_local);
            const double vavg = pgc->globalmax(vavg_local);
            const double nvel = normal_x[n]*uavg + normal_y[n]*vavg;

            if(p->mpirank==0)
            {
                pout[n]<<setprecision(12)
                       <<dist[n][q]<<" "
                       <<p->Xout(xpt[n][q],ypt[n][q])<<" "
                       <<p->Yout(xpt[n][q],ypt[n][q])<<" "
                       <<uavg<<" "
                       <<vavg<<" "
                       <<nvel<<endl;
            }
        }

        pout[n].close();
    }
}

void nhflow_depavg_vel_lineprobe::ini_location(lexer *p, fdm_nhf *d)
{
    for(n=0;n<probenum;++n)
    {
        flag[n] = 1;

        const double xs = p->P144_xs[n];
        const double xe = p->P144_xe[n];
        const double ys = p->P144_ys[n];
        const double ye = p->P144_ye[n];

        const double dx = xe - xs;
        const double dy = ye - ys;

        line_length[n] = sqrt(dx*dx + dy*dy);

        if(line_length[n]>1.0e-14)
        {
            normal_x[n] = -dy/line_length[n];
            normal_y[n] = dx/line_length[n];
        }
        else
        {
            normal_x[n] = 0.0;
            normal_y[n] = 0.0;
        }

        for(q=0;q<point_count[n];++q)
        {
            const double fac = (point_count[n]==1)?0.0:double(q)/double(point_count[n]-1);

            xpt[n][q] = xs + fac*dx;
            ypt[n][q] = ys + fac*dy;
            dist[n][q] = fac*line_length[n];

            iloc[n][q] = p->posc_i(xpt[n][q]);

            if(p->j_dir==0)
            {
                jloc[n][q] = 0;
                j = 0;
                ypt[n][q] = 0.5*p->YP[JP];
            }
            else
            {
                jloc[n][q] = p->posc_j(ypt[n][q]);
            }

            if(!(iloc[n][q]>=0 && iloc[n][q]<p->knox && jloc[n][q]>=0 && jloc[n][q]<p->knoy))
            {
                iloc[n][q] = -1;
                jloc[n][q] = -1;
            }
        }
    }
}
