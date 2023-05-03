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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::distbeach_ini(lexer *p)
{
    double dx,dy,l1,l2,n1x,n1y,n2x,n2y;
    double dist;
    
    for(n=0;n<p->B107;++n)
    {
    Bs[n][0] = p->B107_xs[n];
    Bs[n][1] = p->B107_ys[n];
    
    Be[n][0] = p->B107_xe[n];
    Be[n][1] = p->B107_ye[n];
    }

    for(n=0;n<p->B107;++n)
    {
    dx = Be[n][0]-Bs[n][0]; 
    dy = Be[n][1]-Bs[n][1]; 
    
    n1x = -dy;
    n1y =  dx;
    
    n2x =  dy;
    n2y = -dx;
    
    l1 = sqrt(pow(n1x,2.0) + pow(n1y,2.0));
    l2 = sqrt(pow(n2x,2.0) + pow(n2y,2.0));
    
    l1 = fabs(l1)>1.0e-20?l1:1.0e20;
    l2 = fabs(l2)>1.0e-20?l2:1.0e20;
    
    dist=p->B107_d[n];
    
    //if(p->mpirank==0)
    //cout<<"dist: "<<dist<<" "<<l1<<" "<<l2<<endl;
    
    B1[n][0] = Bs[n][0] + (n1x/l1)*dist; 
    B1[n][1] = Bs[n][1] + (n1y/l1)*dist; 
    
    B2[n][0] = Bs[n][0] + (n2x/l2)*dist; 
    B2[n][1] = Bs[n][1] + (n2y/l2)*dist; 
    
    B3[n][0] = Be[n][0] + (n1x/l1)*dist; 
    B3[n][1] = Be[n][1] + (n1y/l1)*dist; 
    
    B4[n][0] = Be[n][0] + (n2x/l2)*dist; 
    B4[n][1] = Be[n][1] + (n2y/l2)*dist; 
    }
    
    /*
    n=0;
    if(p->mpirank==0)
    cout<<" Beach Polygon  B1:"<<B1[n][0]<<" "<<B1[n][1]<<"  B2: "<<B2[n][0]<<" "<<B2[n][1]<<"  B3: "<<B3[n][0]<<" "<<B3[n][1]<<"  B4: "<<B4[n][0]<<" "<<B4[n][1]<<endl;
    */
}

void iowave::distgen_ini(lexer *p)
{
    double dx,dy,l1,l2,n1x,n1y,n2x,n2y;
    double dist;
    
    
    for(n=0;n<p->B108;++n)
    {
    Gs[n][0] = p->B108_xs[n];
    Gs[n][1] = p->B108_ys[n];
    
    Ge[n][0] = p->B108_xe[n];
    Ge[n][1] = p->B108_ye[n];
    }

    for(n=0;n<p->B108;++n)
    {
    dx = Ge[n][0]-Gs[n][0]; 
    dy = Ge[n][1]-Gs[n][1]; 
    
    n1x = -dy;
    n1y =  dx;
    
    n2x =  dy;
    n2y = -dx;
    
    l1 = sqrt(pow(n1x,2.0) + pow(n1y,2.0));
    l2 = sqrt(pow(n2x,2.0) + pow(n2y,2.0));
    
    l1 = fabs(l1)>1.0e-20?l1:1.0e20;
    l2 = fabs(l2)>1.0e-20?l2:1.0e20;
    
    dist=p->B108_d[n];
    
    //if(p->mpirank==0)
    //cout<<"dist: "<<dist<<" "<<l1<<" "<<l2<<endl;
    
    G1[n][0] = Gs[n][0] + (n1x/l1)*dist; 
    G1[n][1] = Gs[n][1] + (n1y/l1)*dist; 
    
    G2[n][0] = Gs[n][0] + (n2x/l2)*dist; 
    G2[n][1] = Gs[n][1] + (n2y/l2)*dist; 
    
    G3[n][0] = Ge[n][0] + (n1x/l1)*dist; 
    G3[n][1] = Ge[n][1] + (n1y/l1)*dist; 
    
    G4[n][0] = Ge[n][0] + (n2x/l2)*dist; 
    G4[n][1] = Ge[n][1] + (n2y/l2)*dist; 
    }
    
    /*
    n=0;
    if(p->mpirank==0)
    cout<<" Gen Polygon  G1:"<<G1[n][0]<<" "<<G1[n][1]<<"  G2: "<<G2[n][0]<<" "<<G2[n][1]<<"  G3: "<<G3[n][0]<<" "<<G3[n][1]<<"  G4: "<<G4[n][0]<<" "<<G4[n][1]<<endl;
    */
}

int iowave::intriangle(lexer *p, double Ax, double Ay, double Bx, double By, double Cx, double Cy, double x0,double y0)
{
	double Px,Py,Pz;
	double Qx,Qy,Qz;
	double Rx,Ry;
	double PQx,PQy,PQz;
	
	double Mx,My;
	double u,v,w;
					
		Px = x0;
		Py = y0;
		Pz = -10.0;
		
		Qx = x0;
		Qy = y0;
		Qz = +10.0;
		
		
		PQx = Qx-Px;
		PQy = Qy-Py;
		PQz = Qz-Pz; 

		
		// uvw
		Mx = PQy*Pz - PQz*Py;
		My = PQz*Px - PQx*Pz;

		
		u = PQz*(Cx*By - Cy*Bx)
		  + Mx*(Cx-Bx) + My*(Cy-By);
		  
		v = PQz*(Ax*Cy - Ay*Cx)
		  + Mx*(Ax-Cx) + My*(Ay-Cy);
		  
		w = PQz*(Bx*Ay - By*Ax)
		  + Mx*(Bx-Ax) + My*(By-Ay);
		
		
		int check=1;
		if(u==0.0 && v==0.0 && w==0.0)
		check = 0;
		
        if((u>=0.0 && v>=0.0 && w>=0.0) || (u<0.0 && v<0.0 && w<0.0) && check==1)
		return 1;
        
        else
        return 0;
    
}
