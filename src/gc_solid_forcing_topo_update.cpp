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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"mgc1.h"
#include"mgc2.h"
#include"mgc3.h"
#include"mgc4.h"
#include"mgc4a.h"

void ghostcell::solid_forcing_topo_update(lexer *p, fdm *a)
{
    int **cellmem1, **cellmem2, **cellmem3, **cellmem4;
    int cellcount1,cellcount2,cellcount3,cellcount4;
    cellcount1=cellcount2=cellcount3=cellcount4=0;
    int cellmemsize=MAX(p->cellnum,3*p->gcb4_count);
	
	p->Iarray(cellmem1,cellmemsize,4);
	p->Iarray(cellmem2,cellmemsize,4);
	p->Iarray(cellmem3,cellmemsize,4);
	p->Iarray(cellmem4,cellmemsize,4);



//-------------------

    velcell_update(p,a,cellmem1,cellcount1,0.0,0.0,0.5*p->DXM,1);
    velcell_update(p,a,cellmem2,cellcount2,0.0,0.5*p->DXM,0.0,2);
    velcell_update(p,a,cellmem3,cellcount3,0.5*p->DXM,0.0,0.0,3);
	gctopo_pressureupdate(p,a,cellmem4,cellcount4,a->press);

    
    //gcb_velflagio(p,a);

    start1(p,a->u,10);
    start2(p,a->v,11);
    start3(p,a->w,12);
	
	int gcval_press; 
	
    gcval_press=40;  
	
	start4(p,a->press,gcval_press);

    p->del_Iarray(cellmem1,cellmemsize,4);
	p->del_Iarray(cellmem2,cellmemsize,4);
	p->del_Iarray(cellmem3,cellmemsize,4);
	p->del_Iarray(cellmem4,cellmemsize,4);


	count=0;
	LOOP
	++count;

	count=globalisum(count);
    
	
	if(p->mpirank==0)
	cout<<"Topo: active number of cells: "<<count<<endl;
 
    //cout<<p->mpirank<<" p->gcb4_count_topo: "<<p->gcb4_count<<endl;
}
