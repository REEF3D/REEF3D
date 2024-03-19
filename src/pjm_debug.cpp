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

#include"pjm.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
 
void pjm::debug(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w, double alpha)
{/*
    ULOOP
    u(i,j,k) = 0.0;
    
    VLOOP
    v(i,j,k) = 0.0;
    
    WLOOP
    w(i,j,k) = 0.0;
    
    pgc->start1(p,u,999);
	pgc->start2(p,v,999);
	pgc->start3(p,w,999);
    
    
    pip=p->Y50;
    if(p->mpirank==2 || p->mpirank==3)
    LOOP
    {

    
    if(p->flag1[Im1JK]<0 && fabs(u(i-1,j,k))==0.0)
    cout<<p->mpirank<<" U: "<<u(i-1,j,k)<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
    
    if(p->flag2[IJm1K]<0 && fabs(v(i,j-1,k))==0.0)
    cout<<p->mpirank<<" V: "<<v(i,j-1,k)<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
    
    //if(p->flag3[IJKm1]<0 && fabs(w(i,j,k-1))==0.0)
    //cout<<p->mpirank<<" W: "<<w(i,j,k-1)<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
    }
    pip=0;
    
    */
    
    /*
    LOOP
    a->press(i,j,k) = 0.0;
    
    pgc->start4(p,a->press,999);
    
    
    if(p->mpirank==2 || p->mpirank==3)
    {

    ULOOP
    if(fabs(a->press(i,j,k))>0.0)
    cout<<p->mpirank<<" Px: "<<a->press(i,j,k)<<" i: "<<i<<" j: "<<j<<" k: "<<k<<endl;
    
    VLOOP
    if(fabs(a->press(i,j,k))>0.0)
    cout<<p->mpirank<<" Py: "<<a->press(i,j,k)<<" i: "<<i<<" j: "<<j<<" k: "<<k<<endl;
    
    WLOOP
    if(fabs(a->press(i,j,k))>0.0)
    cout<<p->mpirank<<" Pz: "<<a->press(i,j,k)<<" i: "<<i<<" j: "<<j<<" k: "<<k<<endl;
    }*/
    
    pgc->start4(p,a->press,40);
    
    LOOP
    {
        if(p->flag4[Im1JK]==SOLID)
		{
        if(a->press(i-1,j,k)!=a->press(i,j,k))
        cout<<p->mpirank<<" P_i: "<<a->press(i,j,k)<<" P_i-1: "<<a->press(i-1,j,k)<<" P_i-2: "<<a->press(i-2,j,k)<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
		}
		
		if(p->flag4[Ip1JK]==SOLID)
		{
		if(a->press(i+1,j,k)!=a->press(i,j,k))
        {
        cout<<p->mpirank<<" P_i: "<<a->press(i,j,k)<<" P_i+1: "<<a->press(i+1,j,k)<<" P_i+2: "<<a->press(i+2,j,k)
        <<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z();
        cout<<" gcorig_i+1: "<<p->gcorig4[p->mgc4[Ip1JK]-10][3][1]<<" gcorig_i+2: "<<p->gcorig4[p->mgc4[Ip2JK]-10][3][2];
        cout<<" | "<<" | "<<p->flag4[Im3JK]<<" "<<p->flag4[Im2JK]<<" "<<p->flag4[Im1JK]<<" |"<<p->flag4[IJK]<<"| "<<" "
        <<p->flag4[Ip1JK]<<" "<<p->flag4[Ip2JK]<<" "<<p->flag4[Ip3JK]<<endl;
        }
		}
		
		if(p->flag4[IJm1K]==SOLID)
		{
		if(a->press(i,j-1,k)!=a->press(i,j,k))
        cout<<p->mpirank<<" P_j: "<<a->press(i,j,k)<<" P_j-1: "<<a->press(i,j-1,k)<<" P_j-2: "<<a->press(i,j-2,k)<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
		}
		
		if(p->flag4[IJp1K]==SOLID)
		{
		if(a->press(i,j+1,k)!=a->press(i,j,k))
        cout<<p->mpirank<<" P_j: "<<a->press(i,j,k)<<" P_j+1: "<<a->press(i,j+1,k)<<" P_j+2: "<<a->press(i,j+2,k)<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
		}
		
		if(p->flag4[IJKm1]==SOLID)
		{
		if(a->press(i,j,k-1)!=a->press(i,j,k))
        cout<<p->mpirank<<" P_k: "<<a->press(i,j,k)<<" P_k-1: "<<a->press(i,j,k-1)<<" P_k-2: "<<a->press(i,j,k-2)<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
		}
		
		if(p->flag4[IJKp1]==SOLID)
		{
		if(a->press(i,j,k+1)!=a->press(i,j,k))
        cout<<p->mpirank<<" P_k: "<<a->press(i,j,k)<<" P_k+1: "<<a->press(i,j,k+1)<<" P_k+2: "<<a->press(i,j,k+2)<<" | x: "<<p->pos_x()<<" y: "<<p->pos_y()<<" z: "<<p->pos_z()<<endl;
		}
    }

    
    
}
 


