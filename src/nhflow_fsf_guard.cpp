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
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_f::fsf_guard(lexer* p, fdm_nhf* d, ghostcell* pgc, slice& WL, slice &K)
{   
    if(p->A580==1)
    SLICELOOP4
    if(p->flagfsf[IJ]==0)
    {
    
    WL(i,j) = d->depth(i,j);
    K(i,j) = 0.0;
        
    }

    SLICELOOP4
    if(p->flagfsf[IJ]==1)
    {  
        if(p->flagfsf[Im1J]==0)
        WL(i-1,j) = WL(i,j);

        if(p->flagfsf[Ip1J]==0)
        WL(i+1,j) = WL(i,j);
        
        if(p->flagfsf[IJm1]==0)
        WL(i,j-1) = WL(i,j);
        
        if(p->flagfsf[IJp1]==0)
        WL(i,j+1) = WL(i,j);
    }
    
}