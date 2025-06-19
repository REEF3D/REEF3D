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
#include"ioflow.h"

void nhflow_fsf_f::depth_update(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow)
{
    SLICELOOP4
	d->depth(i,j) = p->wd - d->bed(i,j);
    
    SLICELOOP4
    d->WL(i,j) = d->eta(i,j) + d->depth(i,j);
    
    pgc->gcsl_start4(p,d->depth,50);
    pgc->gcsl_start4(p,d->WL,gcval_eta);   
    
    wetdry(p,d,pgc,d->U,d->V,d->W,d->WL); 
}