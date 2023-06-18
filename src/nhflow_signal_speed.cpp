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
FITNESS FOR A PARTICULAR PURPOSE. d->Se[IJK]e the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_signal_speed.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_nhf.h"
#include"patchBC_interface.h"

nhflow_signal_speed::nhflow_signal_speed(lexer* p) 
{

}

nhflow_signal_speed::~nhflow_signal_speed()
{
}

void nhflow_signal_speed::signal_speed_update(lexer* p, ghostcell *pgc, fdm_nhf *d, 
                                        double *Us, double *Un, double *Ve, double *Vw, 
                                        slice &Ds,slice &Dn, slice &De, slice &Dw)
{
    // signal speed x-dir
    ULOOP
    {
    USx = 0.5*(Us[IJK]+Un[IJK]) + sqrt(9.81*Ds(i,j)) - sqrt(9.81*Dn(i,j));
    DSx = 0.5*(sqrt(9.81*Ds(i,j)) + sqrt(9.81*Dn(i,j))) + 0.25*(Us[IJK] - Un[IJK]);
    
    d->Ss[IJK] = MIN(Us[IJK] - sqrt(9.81*Ds(i,j)), USx - DSx);
    d->Sn[IJK] = MAX(Un[IJK] + sqrt(9.81*Dn(i,j)), USx + DSx);
    
    if(p->wet[Ip1J]==0)
    {
    d->Ss[IJK] = Us[IJK] - sqrt(9.81*Ds(i,j));
    d->Sn[IJK] = Us[IJK] + 2.0*sqrt(9.81*Ds(i,j));
    }
    
    if(p->wet[Im1J]==0)
    {
    d->Ss[IJK] = Un[IJK] - 2.0*sqrt(9.81*Dn(i,j));
    d->Sn[IJK] = Un[IJK] + sqrt(9.81*Dn(i,j));
    }
    
    if(p->wet[IJ]==0)
    {
    d->Ss[IJK]=d->Sn[IJK]=0.0;
    USx=0.0;
    }
    }
    
    // signal speed y-dir
    VLOOP
    {
    USy = 0.5*(Ve[IJK]+Vw[IJK]) + sqrt(9.81*De(i,j)) - sqrt(9.81*Dw(i,j));
    DSy = 0.5*(sqrt(9.81*De(i,j)) + sqrt(9.81*Dw(i,j))) + 0.25*(Ve[IJK] - Vw[IJK]);
    
    d->Se[IJK] = MIN(Ve[IJK] - sqrt(9.81*De(i,j)), USy - DSy);
    d->Sw[IJK] = MAX(Vw[IJK] + sqrt(9.81*Dw(i,j)), USy + DSy);
    
    if(p->wet[IJp1]==0)
    {
    d->Se[IJK] = Ve[IJK] - sqrt(9.81*De(i,j));
    d->Sw[IJK] = Ve[IJK] + 2.0*sqrt(9.81*De(i,j));
    }
    
    if(p->wet[IJm1]==0)
    {
    d->Se[IJK] = Vw[IJK] - 2.0*sqrt(9.81*Dw(i,j));
    d->Sw[IJK] = Vw[IJK] + sqrt(9.81*Dw(i,j));
    }
    
    if(p->wet[IJ]==0)
    {
    d->Se[IJK]=d->Sw[IJK]=0.0;
    USy=0.0;
    }
    }
    
}



