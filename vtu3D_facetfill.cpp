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

#include"vtu3D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"cart4.h"

void vtu3D::ggcfacet_fill(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
	int *ggcfacet,**ggcside;
	const double eps = 1.0e-15;
	double nx,ny,nz,norm;
	double vx,vy,vz;
	int q;
	fieldint5 g(p);
	
	LOOP
	g(i,j,k)=0;
	
	p->Iarray(ggcfacet,cart4::ggccount);
	p->Iarray(ggcside,cart4::ggccount,3);
	
	// facet vs ggc
	for(n=0;n<p->facetnum;++n)
	{
    i=p->facet[n][0];
    j=p->facet[n][1];
    k=p->facet[n][2];
	 
	g(i,j,k)=1;
	}
	
	for(n=0;n<cart4::ggccount;n++)
	ggcfacet[n]=0;
	
	for(n=0;n<cart4::ggccount;n++)
	{
	i=cart4::ggc[n][0];
	j=cart4::ggc[n][1];
	k=cart4::ggc[n][2];
	
	if(g(i,j,k)==1)
	ggcfacet[n]=1;
	}
	
	for(n=0;n<cart4::ggccount;n++)
	for(q=0;q<3;++q)
	ggcside[n][q]=0;
	
	
	// x-side
	for(n=0;n<p->facetnum;++n)
	{
    i=p->facet[n][0];
    j=p->facet[n][1];
    k=p->facet[n][2];
	 
	if(p->gcn[n][0]>eps)
	g(i,j,k)=1;
	
	if(p->gcn[n][0]<-eps)
	g(i,j,k)=-1;
	}
	
	for(n=0;n<cart4::ggccount;n++)
	{
	i=cart4::ggc[n][0];
	j=cart4::ggc[n][1];
	k=cart4::ggc[n][2];
	
	if(fabs(g(i,j,k))>0)
	ggcside[n][0]=g(i,j,k);
	}
	
	// y-side
	for(n=0;n<p->facetnum;++n)
	{
    i=p->facet[n][0];
    j=p->facet[n][1];
    k=p->facet[n][2];
	 
	if(p->gcn[n][1]>eps)
	g(i,j,k)=1;
	
	if(p->gcn[n][1]<-eps)
	g(i,j,k)=-1;
	}
	
	for(n=0;n<cart4::ggccount;n++)
	{
	i=cart4::ggc[n][0];
	j=cart4::ggc[n][1];
	k=cart4::ggc[n][2];
	
	if(fabs(g(i,j,k))>0)
	ggcside[n][1]=g(i,j,k);
	}
	
	// z-side
	for(n=0;n<p->facetnum;++n)
	{
    i=p->facet[n][0];
    j=p->facet[n][1];
    k=p->facet[n][2];
	 
	if(p->gcn[n][2]>eps)
	g(i,j,k)=1;
	
	if(p->gcn[n][2]<-eps)
	g(i,j,k)=-1;
	}
	
	for(n=0;n<cart4::ggccount;n++)
	{
	i=cart4::ggc[n][0];
	j=cart4::ggc[n][1];
	k=cart4::ggc[n][2];
	
	if(fabs(g(i,j,k))>0)
	ggcside[n][2]=g(i,j,k);
	}
	
	for(n=0;n<cart4::ggccount;n++)
	if(ggcfacet[n]>0)
	{	
		 nx=ny=nz=norm=0.0;
		 
		 i=cart4::ggc[n][0];
		 j=cart4::ggc[n][1];
		 k=cart4::ggc[n][2];
		 q = 0;// cart4::ggc[n][10]; 

		 if(cart4::ggc[n][4]==1 && ggcside[n][0]==1)
		 {
		 nx = (f(i+1+q,j,k)-f(i+q,j,k))/p->dx;
		 vx = f(i+q,j,k);
		 //if(fabs(nx)>1.000001)
		 //cout<<"1 GGCFACETFILL "<<f(i+1+q,j,k)<<" "<<f(i+q,j,k)<<" nx: "<<nx<<" vx: "<<vx<<" q: "<<q<<endl;
		 }
		 
		 if(cart4::ggc[n][5]==1 && ggcside[n][1]==-q)
		 {
		 ny = (f(i,j-q,k)-f(i,j-1-q,k))/p->dx;
		 vy = f(i,j-q,k);
		 //if(fabs(ny)>1.000001)
		 //cout<<"2 GGCFACETFILL "<<f(i,j-q,k)<<" "<<f(i,j-1-q,k)<<" ny: "<<ny<<" vy: "<<vy<<" q: "<<q<<endl;
		 }

		 if(cart4::ggc[n][6]==1 && ggcside[n][1]==1)
		 {
		 ny = (f(i,j+1+q,k)-f(i,j+q,k))/p->dx;
		 vy = f(i,j+q,k);
		 //if(fabs(ny)>1.000001)
		 //cout<<"3 GGCFACETFILL "<<f(i,j+1+q,k)<<" "<<f(i,j+1+q,k)<<" ny: "<<ny<<" vy: "<<vy<<" q: "<<q<<endl;
		 }

		 if(cart4::ggc[n][7]==1 && ggcside[n][0]==-q)
		 {
		 nx = (f(i-q,j,k)-f(i-1-q,j,k))/p->dx;
		 vx = f(i-q,j,k);
		 //if(fabs(nx)>1.000001)
		 //cout<<"4 GGCFACETFILL "<<f(i-q,j,k)<<" "<<f(i-1-q,j,k)<<" nx: "<<nx<<" vx: "<<vx<<" q: "<<q<<endl;
		 }

		 if(cart4::ggc[n][8]==1 && ggcside[n][2]==1)
		 {
		 nz = (f(i,j,k+1+q)-f(i,j,k+q))/p->dx;
		 vz = f(i,j,k+q);
		 //if(fabs(nz)>1.000001)
		 //cout<<"5 GGCFACETFILL "<<f(i,j,k+1+q)<<" "<<f(i,j,k+q)<<" nz: "<<nz<<" vz: "<<vz<<" q: "<<q<<endl;
		 }

		 if(cart4::ggc[n][9]==1 && ggcside[n][2]==-q)
		 {
		 nz = (f(i,j,k-q)-f(i,j,k-1-q))/p->dx;
		 vz = f(i,j,k-q);
		 //if(fabs(nz)>1.000001)
		 //cout<<"6 GGCFACETFILL "<<f(i,j,k-q)<<" "<<f(i,j,k-1-q)<<" nz: "<<nz<<" vz: "<<vz<<" q: "<<q<<endl;
		 }
		 
		 norm = sqrt(nx*nx + ny*ny + nz*nz);
		 
		 //if(fabs(nx)>1.000001||fabs(ny)>1.000001||fabs(nz)>1.000001)
		 //cout<<"# GGCFACETFILL "<<nx<<"  "<<ny<<" "<<nz<<"  "<<norm<<" | "<<vx<<"  "<<vy<<" "<<vz<<"  "<<endl<<endl;
		 
		 nx/=norm>1.0e-20?norm:1.0e20;
		 ny/=norm>1.0e-20?norm:1.0e20;
		 nz/=norm>1.0e-20?norm:1.0e20;
		 
		 pip=4;
		 f(i,j,k) = vx*(1.0-fabs(nx)) + vy*(1.0-fabs(ny)) + vz*(1.0-fabs(nz));
		 pip=0;
		 
		 //cout<<"GGCFACETFILL "<<nx<<"  "<<ny<<" "<<nz<<"  "<<norm<<" | "<<vx<<"  "<<vy<<" "<<vz<<"  "<<endl;
		 
	 }
	
	p->del_Iarray(ggcfacet,cart4::ggccount);
	p->del_Iarray(ggcside,cart4::ggccount,3);
}