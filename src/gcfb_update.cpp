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

void ghostcell::gcfb_update(lexer *p, fdm *a)
{
    if(p->mpirank==0)
    cout<<"6DOF: update grid..."<<endl;

    int **cellmem1, **cellmem2, **cellmem3, **cellmem4;
    int cellcount1,cellcount2,cellcount3,cellcount4;
    cellcount1=cellcount2=cellcount3=cellcount4=0;
    int cellmemsize=p->cellnum;

	mgc1 m1(p);
	mgc2 m2(p);
	mgc3 m3(p);
	mgc4 m4(p);

	p->Iarray(cellmem1,cellmemsize,8);
	p->Iarray(cellmem2,cellmemsize,8);
	p->Iarray(cellmem3,cellmemsize,8);
	p->Iarray(cellmem4,cellmemsize,8);

    gcfb_buildflag(p,a,cellmem4,cellcount4);

    flagx(p,p->flag4);
    gcb_remove(p,a);
    gcfb_seed(p,a);
    gcfb_dist(p,a);

    gcfb_velflag1(p,a,cellmem1,cellcount1);
    gcfb_velflag2(p,a,cellmem2,cellcount2);
    gcfb_velflag3(p,a,cellmem3,cellcount3);
	
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
    
//
    sizeM_update(p,a);
  
    m1.make_dgc(p);
    m2.make_dgc(p);
    m3.make_dgc(p);
    m4.make_dgc(p);
    
    m1.fill_dgc(p);
    m2.fill_dgc(p);
    m3.fill_dgc(p);
    m4.fill_dgc(p);
  
    column_pt_resize(p,a);
    column_pt4_update(p,a); 

//-------------------

    gcfb_velupdate(p,a,cellmem1,cellcount1,0.0,0.0,0.5*p->DXM,1);
    gcfb_velupdate(p,a,cellmem2,cellcount2,0.0,0.5*p->DXM,0.0,2);
    gcfb_velupdate(p,a,cellmem3,cellcount3,0.5*p->DXM,0.0,0.0,3);
	//gcfb_scalarupdate(p,a,cellmem4,cellcount4,a->press);
	
//-------------------
	

    start1(p,a->u,10);
    start2(p,a->v,11);
    start3(p,a->w,12);
	
	int gcval_press; 

    gcval_press=40;  

    start4(p,a->press,gcval_press);

    p->del_Iarray(cellmem1,cellmemsize,8);
	p->del_Iarray(cellmem2,cellmemsize,8);
	p->del_Iarray(cellmem3,cellmemsize,8);
	p->del_Iarray(cellmem4,cellmemsize,8);
	
	
    //gcparax_test(a,p,4);
    
	count=0;
	LOOP
	++count;
	
	count=globalisum(count);
	
}
