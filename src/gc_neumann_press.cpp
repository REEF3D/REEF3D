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

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include"cpt.h"

void ghostcell::neumann_press(lexer *p,field& f,double dist,int gcv, int bc, int cs)
{
	double ui,vi,wi,un,vn,wn;
    
    double dx;
    
    if(cs==1||cs==4)
    dx = p->DXP[IP];
    
    if(cs==2||cs==3)
    dx = p->DYP[JP];
    
    if(cs==5||cs==6)
    dx = p->DZP[KP];
	
	wallvalue=0.0;
	weight=1.0;
	margin=3;
		
	if(bc==41 || bc==42 || bc==43)
	{
		
	ui = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
	vi = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
	wi = p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi;
	
	un = p->ufbn + (p->pos_z()-p->zgn)*p->qfbn - (p->pos_y()-p->ygn)*p->rfbn;
	vn = p->vfbn + (p->pos_x()-p->xgn)*p->rfbn - (p->pos_z()-p->zgn)*p->pfbn;
	wn = p->wfbn + (p->pos_y()-p->ygn)*p->pfbn - (p->pos_x()-p->xgn)*p->qfbn;

	if(cs==1)
	wallvalue =  (ui-un)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==2)
	wallvalue = -  (vi-vn)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==3)
	for(q=0;q<margin;++q)
	wallvalue =  (vi-vn)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==4)
	for(q=0;q<margin;++q)
	wallvalue = -  (ui-un)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==5)
	for(q=0;q<margin;++q)
	wallvalue =  (wi-wn)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==6)
	for(q=0;q<margin;++q)
	wallvalue = -  (wi-wn)/(p->dt>1.0e-20?p->dt:1.0e20);
	
	wallvalue *= a->ro(i,j,k);
	}	
	
	ys=1;
    if(dist>0.5*dx*(1.0-1.0e-6) && dist<0.5*dx*(1.0+1.0e-6))
    ys=0;
	
    if(ys==1)
    {

//fill pos[]
	for(m=0;m<=orderpress-3;m++)
	pos[m]=-dx*double(orderpress-m-2)-0.5*dx;

	pos[orderpress-2]=-0.5*dx;
	pos[orderpress-1]=dist;

	for(m=0;m<margin;m++)
	x[m] = dx*double(m+2-ys)-0.5*dx;

//fill y[]
	if(cs==1 )
	for(m=0;m<orderpress-1;m++)
	y[m]=(f(i+orderpress-m-1,j,k)-f(i+orderpress-m-2,j,k))/dx;

	if(cs==2)
	for(m=0;m<orderpress-1;m++)
	y[m]=(f(i,j-orderpress+m+2,k)-f(i,j-orderpress+m+1,k))/dx;

	if(cs==3)
	for(m=0;m<orderpress-1;m++)
	y[m]=(f(i,j+orderpress-m-1,k)-f(i,j+orderpress-m-2,k))/dx;

	if(cs==4)
	for(m=0;m<orderpress-1;m++)
	y[m]=(f(i-orderpress+m+2,j,k)-f(i-orderpress+m+1,j,k))/dx;

	if(cs==5)
	for(m=0;m<orderpress-1;m++)
	y[m]=(f(i,j,k+orderpress-m-1)-f(i,j,k+orderpress-m-2))/dx;

	if(cs==6)
	for(m=0;m<orderpress-1;m++)
	y[m]=(f(i,j,k-orderpress+m+2)-f(i,j,k-orderpress+m+1))/dx;

	y[orderpress-1]=wallvalue;

    if(ys==1 && dist<gamma*dx)
    {
    imagepoint(p,f,x_ip,val_ip,dist,cs);

    pos[orderpress-2] = x_ip;
    y[orderpress-2] = val_ip;
    }

	for(q=0; q<margin; ++q)
	{
	    y[orderpress+q]=0.0;

		for(m=0;m<orderpress;m++)
		{
			weight=1.0;
			for(n=0;n<orderpress;++n)
			{
			if(m!=n)
			weight*=(x[q]-pos[n])/(pos[m]-pos[n]+1.0e-20);
			}
		y[orderpress+q]+=weight*y[m];
		}
	}
	
	/*
	if(cs==5)
	{
	cout<<dist<<" . "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<" || "<<y[0]<<" "<<y[1]<<" "<<y[2]<<" "<<y[3]<<" "<<y[4]<<" "<<y[5]<<" || ";
	cout<<f(i,j,k)<<" . ";
	for(q=0;q<margin;++q)
	cout<<f(i,j,k) - double(q+1)*dx*y[orderpress+q]<<" ";
	
	cout<<endl;
	}*/

// write extrapolated ghost cell values into f()

	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=f(i,j,k) - (double(q)*dx + dx)*y[orderpress+q-1+ys];

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=f(i,j,k) + (double(q)*dx + dx)*y[orderpress+q-1+ys];

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=f(i,j,k) - (double(q)*dx + dx)*y[orderpress+q-1+ys];

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=f(i,j,k) + (double(q)*dx + dx)*y[orderpress+q-1+ys];

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=f(i,j,k) - (double(q)*dx + dx)*y[orderpress+q-1+ys];

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=f(i,j,k) + (double(q)*dx + dx)*y[orderpress+q-1+ys];

    }


    if(ys==0)
    {
    if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=f(i,j,k);

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=f(i,j,k);

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=f(i,j,k);

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=f(i,j,k);

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=f(i,j,k) - double(q+1)*dx*wallvalue;

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=f(i,j,k) + double(q+1)*dx*wallvalue;
    }
}

