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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include "grid.h"
#include "ghostcell.h"
#include "resize.h"
#include <cmath>
#include <vector>
#include <algorithm>

void grid::assign_margin()
{
    imax=knox+2*margin;
    jmax=knoy+2*margin;
    kmax=knoz+2*margin;
    kmaxF=knoz+1+2*margin;

    imin=-margin;
    jmin=-margin;
    kmin=-margin;
    
    //global_orig_x = global_orig_y = 0.0;
    //alpha_grid = 0.0;
}

void grid::sigma_coord_ini()
{
    double L, ZN0temp;

    L = ZN[knoz+marge] - ZN[0+marge];

    ZN0temp = ZN[0+marge];

    for(k=-marge;k<knoz+marge;++k)
    ZN[KP] = (ZN[KP]-ZN0temp)/L;
}

void grid::gridspacing(ghostcell *pgc)
{
    resize_class resizer;

    resizer.Darray(XP,knox+1+4*marge);
    resizer.Darray(YP,knoy+1+4*marge);
    resizer.Darray(ZP,knoz+1+4*marge);

    resizer.Darray(DXN,knox+1+4*marge);
    resizer.Darray(DYN,knoy+1+4*marge);
    resizer.Darray(DZN,knoz+1+4*marge);

    resizer.Darray(DXP,knox+1+4*marge);
    resizer.Darray(DYP,knoy+1+4*marge);
    resizer.Darray(DZP,knoz+1+4*marge);

    resizer.Darray(ZSN,imax*jmax*(kmax+1));
    resizer.Darray(ZSP,imax*jmax*kmax);

    // XP,YP,ZP
    for(i=-marge;i<knox+marge;++i)
    XP[IP] = 0.5*(XN[IP]+XN[IP1]);

    for(j=-marge;j<knoy+marge;++j)
    YP[JP] = 0.5*(YN[JP]+YN[JP1]);

    for(k=-marge;k<knoz+marge;++k)
    ZP[KP] = 0.5*(ZN[KP]+ZN[KP1]);

    //dx
    for(i=-marge;i<knox+marge;++i)
    DXN[IP] = XN[IP1]-XN[IP];

    for(j=-marge;j<knoy+marge;++j)
    DYN[JP] = YN[JP1]-YN[JP];

    for(k=-marge;k<knoz+marge;++k)
    DZN[KP] = ZN[KP1]-ZN[KP];

    // dxn

    for(i=-marge;i<knox+marge;++i)
    DXP[IP] = 0.5*(XN[IP2]+XN[IP1]) - 0.5*(XN[IP1]+XN[IP]);

    for(j=-marge;j<knoy+marge;++j)
    DYP[JP] = 0.5*(YN[JP2]+YN[JP1]) - 0.5*(YN[JP1]+YN[JP]);

    for(k=-marge;k<knoz+marge;++k)
    DZP[KP] = 0.5*(ZN[KP2]+ZN[KP1]) - 0.5*(ZN[KP1]+ZN[KP]);

    // every cell spacing once and summed in the same order whatever the partition: the x spacings
    // from the ranks at the lower y and z end of the domain, the y spacings from those at the lower
    // x and z end, the z spacings from those at the lower x and y end, gathered and summed in the
    // order of their coordinates (a rank holds all cells of the directions that are not split;
    // summed on every rank, the mean depended on the partition, and partial sums on the rounding)
    i=j=k=0;
    const double x0 = pgc->globalmin(XN[IP]);
    const double y0 = pgc->globalmin(YN[JP]);
    const double z0 = pgc->globalmin(ZN[KP]);
    const bool lowx = XN[IP]-x0 <= 1.0e-9*MAX(fabs(x0),1.0);
    const bool lowy = YN[JP]-y0 <= 1.0e-9*MAX(fabs(y0),1.0);
    const bool lowz = ZN[KP]-z0 <= 1.0e-9*MAX(fabs(z0),1.0);

    // (direction, coordinate, spacing) of the cells this rank contributes
    std::vector<double> loc;
    if(lowy && lowz)
    for(i=0;i<knox;++i)
    {
        loc.push_back(0.0); loc.push_back(XP[IP]); loc.push_back(DXP[IP]);
    }

    if(j_dir==1 && lowx && lowz)
    for(j=0;j<knoy;++j)
    {
        loc.push_back(1.0); loc.push_back(YP[JP]); loc.push_back(DYP[JP]);
    }

    if(lowx && lowy)
    for(k=0;k<knoz;++k)
    {
        loc.push_back(2.0); loc.push_back(ZP[KP]); loc.push_back(DZP[KP]);
    }

    std::vector<double> all;
    const bool lred = pgc->set_local(false);
    pgc->set_local(lred);
    int np=1;
    if(!lred && pgc->mpi_comm!=MPI_COMM_NULL)
    MPI_Comm_size(pgc->mpi_comm,&np);

    if(np>1)
    {
        int n = int(loc.size());
        std::vector<int> cnt(np), dsp(np,0);
        MPI_Allgather(&n,1,MPI_INT,cnt.data(),1,MPI_INT,pgc->mpi_comm);
        for(int q=1; q<np; ++q)
        dsp[q] = dsp[q-1]+cnt[q-1];
        all.resize(dsp[np-1]+cnt[np-1]);
        MPI_Allgatherv(loc.data(),n,MPI_DOUBLE,all.data(),cnt.data(),dsp.data(),MPI_DOUBLE,pgc->mpi_comm);
    }
    else
    all = loc;

    std::vector<int> ord(all.size()/3);
    for(size_t q=0; q<ord.size(); ++q)
    ord[q] = int(q);
    std::stable_sort(ord.begin(),ord.end(),[&](int u, int w)
    {
        if(all[3*u]!=all[3*w])
        return all[3*u]<all[3*w];
        return all[3*u+1]<all[3*w+1];
    });

    DXM = DXD = DYD = 0.0;

    int count=0;
    int xcount=0;
    int ycount=0;

    for(int q : ord)
    {
        const double d = all[3*q+2];
        DXM += d;
        ++count;

        if(all[3*q]==0.0)
        {
        DXD += d;
        ++xcount;
        }

        if(all[3*q]==1.0)
        {
        DYD += d;
        ++ycount;
        }
    }

    DXM /= double(count);
    DXD /= double(xcount);
    DYD /= double(ycount);

    DXM = pgc->globalmin(DXM);
    DXD = pgc->globalmin(DXD);
    DYD = pgc->globalmin(DYD);

    
}
