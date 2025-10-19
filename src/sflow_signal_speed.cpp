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
FITNESS FOR A PARTICULAR PURPOSE. b->Se(i,j)e the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sflow_signal_speed.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"patchBC_interface.h"

sflow_signal_speed::sflow_signal_speed(lexer* p) 
{

}

sflow_signal_speed::~sflow_signal_speed()
{
}

void sflow_signal_speed::signal_speed_update(lexer* p, ghostcell *pgc, fdm2D *b, 
                                        slice &Us, slice &Un, slice &Ve, slice &Vw, 
                                        slice &Ds, slice &Dn, slice &De, slice &Dw)
{
    // signal speed x-dir
    SLICELOOP1
    {
    USx = 0.5*(Us(i,j)+Un(i,j)) + sqrt(9.81*Ds(i,j)) - sqrt(9.81*Dn(i,j));
    DSx = 0.5*(sqrt(9.81*Ds(i,j)) + sqrt(9.81*Dn(i,j))) + 0.25*(Us(i,j) - Un(i,j));
    
    if(p->wet[IJ]==1 && p->wet[Ip1J]==1)
    {
    b->Ss(i,j) = MIN(Us(i,j) - sqrt(9.81*Ds(i,j)), USx - DSx);
    b->Sn(i,j) = MAX(Un(i,j) + sqrt(9.81*Dn(i,j)), USx + DSx);
    b->SSx(i,j) = USx;
    }
    
    else
    if(p->wet[IJ]==0 && p->wet[Ip1J]==1) // left dry
    {
    b->Ss(i,j) = Un(i,j) - 2.0*sqrt(9.81*Dn(i,j));
    b->Sn(i,j) = Un(i,j) + sqrt(9.81*Dn(i,j));
    b->SSx(i,j) = b->Ss(i,j);
    }
    
    else
    if(p->wet[IJ]==1 && p->wet[Ip1J]==0)// right dry
    {
    b->Ss(i,j) = Us(i,j) - sqrt(9.81*Ds(i,j));
    b->Sn(i,j) = Us(i,j) + 2.0*sqrt(9.81*Ds(i,j));
    b->SSx(i,j) = b->Sn(i,j);
    }
    
    else
    if(p->wet[IJ]==0 && p->wet[Ip1J]==0) 
    {
    b->Ss(i,j) = 0.0;
    b->Sn(i,j) = 0.0;
    b->SSx(i,j) = 0.0;
    }
    }
    
    // signal speed y-dir
    if(p->j_dir==1)
    SLICELOOP2
    {
    USy = 0.5*(Ve(i,j)+Vw(i,j)) + sqrt(9.81*De(i,j)) - sqrt(9.81*Dw(i,j));
    DSy = 0.5*(sqrt(9.81*De(i,j)) + sqrt(9.81*Dw(i,j))) + 0.25*(Ve(i,j) - Vw(i,j));
    
    if(p->wet[IJ]==1 && p->wet[IJp1]==1)
    {
    b->Se(i,j) = MIN(Ve(i,j) - sqrt(9.81*De(i,j)), USy - DSy);
    b->Sw(i,j) = MAX(Vw(i,j) + sqrt(9.81*Dw(i,j)), USy + DSy);
    b->SSy(i,j) = USy;
    }

    else
    if(p->wet[IJ]==0 && p->wet[IJp1]==1)
    {
    b->Se(i,j) = Vw(i,j) - 2.0*sqrt(9.81*Dw(i,j));
    b->Sw(i,j) = Vw(i,j) + sqrt(9.81*Dw(i,j));
    b->SSy(i,j) = b->Se(i,j);
    }
    
    else
    if(p->wet[IJ]==1 && p->wet[IJp1]==0)
    {
    b->Se(i,j) = Ve(i,j) - sqrt(9.81*De(i,j));
    b->Sw(i,j) = Ve(i,j) + 2.0*sqrt(9.81*De(i,j));
    b->SSy(i,j) = b->Sw(i,j);
    }
    
    else
    if(p->wet[IJ]==0 && p->wet[IJp1]==0)
    {
    b->Se(i,j) = 0.0;
    b->Sw(i,j) = 0.0;
    b->SSy(i,j) = 0.0;
    }
    
    }
}
