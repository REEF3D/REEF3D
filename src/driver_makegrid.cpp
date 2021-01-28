/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"driver.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mgc1.h"
#include"mgc2.h"
#include"mgc3.h"
#include"mgc4.h"
#include"mgc4a.h"
#include"mgc6.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"cart4a.h"

void driver::makegrid(lexer *p, ghostcell *pgc)
{	
	mgc1 m1(p);
	mgc2 m2(p);
	mgc3 m3(p);
	mgc4 m4(p);
	mgc4a m4a(p);
    mgc6 m6(p);
    
    if(p->mpirank==0)
    {
    i=0;
    j=0;
    k=6;
    
    cout<<i<<" "<<" FLAG PERIODX 1: "<<p->flag4[IJK]<<" "<<p->flag4[Im1JK]<<" "<<p->flag4[Im2JK]<<" "<<p->flag4[Im3JK]<<endl;
    }
    
    if(p->mpirank==5)
    {
    i=p->knox-1;
    j=0;
    k=6;
    
    cout<<i<<" "<<" FLAG1 PERIODX 4: "<<p->flag4[IJK]<<" "<<p->flag4[Ip1JK]<<" "<<p->flag4[Ip2JK]<<" "<<p->flag4[Ip3JK]<<endl;
    }
    
	pgc->flagx(p,p->flag1);
    pgc->flagx(p,p->flag2);
    pgc->flagx(p,p->flag3);
    pgc->flagx(p,p->flag4);
    pgc->flagx(p,p->flag);
	pgc->gcxupdate(p);
    
    
    if(p->mpirank==0)
    {
    i=0;
    j=0;
    k=6;
    
    cout<<i<<" "<<" FLAG PERIODX 1: "<<p->flag4[IJK]<<" "<<p->flag4[Im1JK]<<" "<<p->flag4[Im2JK]<<" "<<p->flag4[Im3JK]<<endl;
    }
    
    if(p->mpirank==5)
    {
    i=p->knox-1;
    j=0;
    k=6;
    
    cout<<i<<" "<<" FLAG1 PERIODX 4: "<<p->flag4[IJK]<<" "<<p->flag4[Ip1JK]<<" "<<p->flag4[Ip2JK]<<" "<<p->flag4[Ip3JK]<<endl;
    }
    
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001"<<endl;
	
	m1.makemgc(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001a"<<endl;
    m1.fillgcb(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001b"<<endl;
    m1.extragcb(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001c"<<endl;
    m1.mgcsetup(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001d"<<endl;
    m1.fillmgc(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001e"<<endl;
    m1.gcdirfill(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  002f"<<endl;
	
    m2.makemgc(p);
    m2.fillgcb(p);
    m2.extragcb(p);
    m2.mgcsetup(p);
    m2.fillmgc(p);
    m2.gcdirfill(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  003"<<endl;

    m3.makemgc(p);
    m3.fillgcb(p);
    m3.extragcb(p);
    m3.mgcsetup(p);
    m3.fillmgc(p);
    m3.gcdirfill(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  004"<<endl;

    m4.makemgc(p);
    m4.mgcsetup(p);
    m4.fillmgc(p);
    m4.gcdirfill(p);
	m4.gcsidefill(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  005"<<endl;

    m4a.makemgc(p);
    m4a.fillgcb(p);
    m4a.mgcsetup(p);
    m4a.fillmgc(p);
    m4a.gcdirfill(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  006"<<endl;
    
    m6.makemgc(p);
    m6.mgcsetup(p);
    m6.fillmgc(p);
    m6.gcdirfill(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  007"<<endl;

	
	m1.make_ggc(p);
    m1.fill_ggc(p);
	m2.make_ggc(p);
    m2.fill_ggc(p);
	m3.make_ggc(p);
    m3.fill_ggc(p);
    m4.make_ggc(p);
    m4.fill_ggc(p);
    m4a.make_ggc(p);
    m4a.fill_ggc(p);
    
    if(p->mpirank==0)
    cout<<"MAKEGRID  001"<<endl;
    
    m1.make_dgc(p);
    m2.make_dgc(p);
    m3.make_dgc(p);
    m4.make_dgc(p);
    
    m1.fill_dgc(p);
    m2.fill_dgc(p);
    m3.fill_dgc(p);
    m4.fill_dgc(p);
    
    p->vecsize(pgc);
}
	
void driver::makegrid_cds()
{	
	pgc->sizeM_update(p,a);
    
    pgc->column_pt4_update(p,a);
    pgc->column_pt4a_update(p,a);
    pgc->column_pt6_update(p,a);
}
	
