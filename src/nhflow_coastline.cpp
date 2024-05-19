/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_coastline.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include"sliceint.h"

nhflow_coastline::nhflow_coastline(lexer* p) :  ddweno_f_nug(p), frk1(p),frk2(p),L(p),dt(p),wet_n(p)
{
    time_preproc(p); 
}

nhflow_coastline::~nhflow_coastline()
{
}

void nhflow_coastline::start(lexer *p, ghostcell *pgc, slice &coastline, int *wet, sliceint &wet_n)
{
    if(p->count==0)
    {
        SLICELOOP4
        {
            if(wet[IJ]==0)
            coastline(i,j)=-1.0;
            
            if(wet[IJ]==1)
            coastline(i,j)=1.0;
   
        }
        reini(p,pgc,coastline);
    }
}




