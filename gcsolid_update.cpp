/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"mgc1.h"
#include"mgc2.h"
#include"mgc3.h"
#include"mgc4.h"
#include"mgc4a.h"

void ghostcell::solid_update(lexer *p, fdm *a)
{
    if(p->mpirank==0)
    cout<<"Solid: update grid..."<<endl;
    
    int cellcount1,cellcount2,cellcount3,cellcount4;
	
	mgc1 m1(p);
	mgc2 m2(p);
	mgc3 m3(p);
	mgc4 m4(p);

    gcsolid_buildflag(p,a,cellcount4);
    
    flagx(p,p->flag4);
    flagx(p,p->flag);
    gcsolid_gcb_remove(p,a);
    gcsolid_gcb_seed(p,a);
    gcsolid_gcb_dist(p,a);


    gcsolid_velflag1(p,a,cellcount1);
    gcsolid_velflag2(p,a,cellcount2);
    gcsolid_velflag3(p,a,cellcount3);

    flagx(p,p->flag1);
    flagx(p,p->flag2);
    flagx(p,p->flag3);
    
    gcxupdate(p);
    
    m1.fillgcb(p);
    m1.extragcb(p);
    m1.mgcsetup(p);
    m1.fillmgc(p);
    m1.gcdirfill(p);

    m2.fillgcb(p);
    m2.extragcb(p);
    m2.mgcsetup(p);
    m2.fillmgc(p);
    m2.gcdirfill(p);

	m3.fillgcb(p);
    m3.extragcb(p);
    m3.mgcsetup(p);
    m3.fillmgc(p);
    m3.gcdirfill(p);

    m4.mgcsetup(p);
    m4.fillmgc(p);
    m4.gcdirfill(p);
	m4.gcsidefill(p);

    m1.fill_ggc(p);
    m2.fill_ggc(p);
    m3.fill_ggc(p);
    m4.fill_ggc(p);
	
	ndflag_update(p);
    

    sizeM_update(p,a);
	
	dgcini1(p);
	dgcini2(p);
	dgcini3(p);
	dgcini4(p);

    column_pt_resize(p,a);
    
    if(p->N40>10 || p->D20>=2)
    {
    column_pt1_update(p,a);
    column_pt2_update(p,a);
    column_pt3_update(p,a);
    }
    column_pt4_update(p,a);

   // a->press.resize(p);


	count=0;
	LOOP
	++count;

	count=globalisum(count);
	
	if(p->mpirank==0)
	cout<<"Solid: active number of cells: "<<count<<endl;
    
    //cout<<p->mpirank<<" GCB4_COUNT:  "<<p->gcb4_count<<endl;
    
    
    if(p->mpirank==0)
    cout<<endl;
}
