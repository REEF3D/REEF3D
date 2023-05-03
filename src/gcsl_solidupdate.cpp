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
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"mgc1.h"
#include"mgc2.h"
#include"mgc3.h"
#include"mgc4.h"
#include"mgc4a.h"

void ghostcell::gcsl_solidupdate(lexer *p)
{
    /*
    if(p->mpirank==0)
    cout<<"Solid Slice: update grid..."<<endl;
    
    int cellcount1,cellcount2,cellcount3,cellcount4;

	mgc4 m4(p);

    //gcsolid_buildflag(p,a,cellcount4);
    
    flagx(p,p->flag4);
    flagx(p,p->flag);
    gcsolid_gcb_remove(p,a);
    gcsolid_gcb_seed(p,a);
    gcsolid_gcb_dist(p,a);
    
    p->gridini_patchBC();	


    gcsolid_velflag1(p,a,cellcount1);
    gcsolid_velflag2(p,a,cellcount2);
    gcsolid_velflag3(p,a,cellcount3);

    flagx(p,p->flag1);
    flagx(p,p->flag2);
    flagx(p,p->flag3);
    
    
    
    
    
        // 2D
    
    gcsl_tpflag(p);    
    gcslflagx(p,p->flagslice4);
    
    mgcslice4 msl4(p);
    
    msl4.makemgc(p);
    msl4.gcb_seed(p);
    msl4.mgcsetup(p);
    msl4.fillmgc(p);
    msl4.gcdirfill(p);
    
    msl4.make_ggc(p);
    msl4.fill_ggc(p);
    
    gcsl_setbc4(p);
    gcsl_setbcio(p);
    
	dgcslini4(p);
    
    
    
    
    // 2D cds
    
    
    sizeS_update(p);
    gcxslupdate(p);
    column2D_pt4_update(p,c->C4); */
}
