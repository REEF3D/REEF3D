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

#include"patch_obj.h"
#include"patchBC_codes.h"
#include"lexer.h"
#include<cmath>

patch_obj::patch_obj(lexer *p, int ID_ini) 
{
    ID = ID_ini;
    kind = PATCH_OUTLET;
    gcb_flag = PATCH_OUTLET;
    
    gcb_count=0;
    gcb=nullptr;
    
    Q_flag=0;
    Q=0.0;
    Uq=0.0;
    
    Uio_flag=0;
    Uio=0.0;
    
    velcomp_flag=0;
    U=V=W=0.0;
    
    angle_flag=0;
    alpha=0.0;
    
    dir_flag=0;
    dirx=diry=dirz=0.0;
    
    pressure_flag=0;
    pressure=0.0;
    
    pio_flag=0;
    
    waterlevel_flag=0;
    waterlevel=0.0;
    
    hydroQ_flag=0;
    hydroQ=nullptr;
    hydroQ_count=0;
    
    hydroFSF_flag=0;
    hydroFSF=nullptr;
    hydroFSF_count=0;
    
    Q0 = 0.0;
    U0 = 0.0;
    A0 = 0.0;
    h0 = 0.0;
}

patch_obj::~patch_obj()
{
}

void patch_obj::velocity(int cs, double Un, double &uval, double &vval, double &wval) const
{
    // B 415: the components as given
    if(velcomp_flag==1)
    {
    uval = U;
    vval = V;
    wval = W;
    return;
    }
    
    double ex=0.0, ey=0.0, ez=0.0;
    
    if(cs==1 || cs==4)
    ex=1.0;
    
    if(cs==2 || cs==3)
    ey=1.0;
    
    if(cs==5 || cs==6)
    ez=1.0;
    
    double dx=ex, dy=ey, dz=ez;
    
    // B 416: horizontal angle, counter-clockwise from the face normal (side faces only)
    if(angle_flag==1)
    {
        if(cs==1 || cs==4)
        {
        dx = cos(alpha);
        dy = sin(alpha);
        dz = 0.0;
        }
        
        if(cs==2 || cs==3)
        {
        dx = -sin(alpha);
        dy = cos(alpha);
        dz = 0.0;
        }
    }
    
    // B 417: direction vector
    if(dir_flag==1)
    {
    dx = dirx;
    dy = diry;
    dz = dirz;
    }
    
    // along the flow direction, scaled so that the component normal to the face is Un
    // (|dir . e| >= 0.1 is checked at the initialisation)
    double dn = dx*ex + dy*ey + dz*ez;
    
    uval = Un*dx/dn;
    vval = Un*dy/dn;
    wval = Un*dz/dn;
}
