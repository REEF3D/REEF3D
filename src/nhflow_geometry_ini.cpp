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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"6DOF.h"
#include"nhflow_reinidisc_fsf.h"

void nhflow_geometry::geometry_ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);
    
    pgc->start5V(p,p->ZSP,1);
    
    
    //DSM
    int num=0;
    DSM=0.0;
    
    if(p->j_dir==0)
    SLICELOOP4
    {
    DSM += p->DXN[IP];
        
    ++num;
    }
    
    if(p->j_dir==1)
    SLICELOOP4
    {
    DSM += 0.5*(p->DXN[IP] + p->DYN[JP]);
        
    ++num;
    }
    
    pgc->globalsum(DSM);
    pgc->globalisum(num);
    
    DSM = DSM/double(num);
}

