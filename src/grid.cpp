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

void grid::assign_margin()
{
    imax=knox+2*margin;
    jmax=knoy+2*margin;
    kmax=knoz+2*margin;
    kmaxF=knoz+1+2*margin;

    imin=-margin;
    jmin=-margin;
    kmin=-margin;
    
    global_orig_x = global_orig_y = 0.0;
    alpha_grid = 0.0;
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

    DXM = DXD = DYD = 0.0;

    int count=0;
    int xcount=0;
    int ycount=0;

    for(i=0;i<knox;++i)
    {
        DXM += DXP[IP];
        DXD += DXP[IP];
        ++count;
        ++xcount;
    }

    if(j_dir==1)
    for(j=0;j<knoy;++j)
    {
        DXM += DYP[JP];
        DYD += DYP[JP];
        ++count;
        ++ycount;
    }

    for(k=0;k<knoz;++k)
    {
        DXM += DZP[KP];
        ++count;
    }

    count = pgc->globalisum(count);
    xcount = pgc->globalisum(xcount);
    ycount = pgc->globalisum(ycount);

    DXM = pgc->globalsum(DXM);
    DXD = pgc->globalsum(DXD);
    DYD = pgc->globalsum(DYD);

    DXM /= double(count);
    DXD /= double(xcount);
    DYD /= double(ycount);

    DXM = pgc->globalmin(DXM);
    DXD = pgc->globalmin(DXD);
    DYD = pgc->globalmin(DYD);

    
}
