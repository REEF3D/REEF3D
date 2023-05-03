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
#include"ghostcell.h"
#include"field.h"
#include"math.h"
#include"density_f.h"

void ghostcell::dirichlet_ortho(lexer *p,field& f,double dist,int gcv, int bc, int cs)
{
	wallvalue=0.0;
	weight=1.0;
	margin=3;
    
    double dx;
    
    if(cs==1||cs==4)
    dx = p->DXP[IP];
    
    if(cs==2||cs==3)
    dx = p->DYP[JP];
    
    if(cs==5||cs==6)
    dx = p->DZP[KP];
	
	ys=1;
    if(dist>p->DXM*(1.0-1.0e-6) && dist<p->DXM*(1.0+1.0e-6))
    ys=0;
    

//fill pos[]
	for(m=0;m<=orderdir-3;m++)
	pos[m]=-dx*double(orderdir-m-2);

	pos[orderdir-2]=0.0;
	pos[orderdir-1]=dist;

	for(m=0;m<margin;m++)
	x[m]=dx*double(m+2-ys);

//fill y[]
	if(cs==1 )
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i+orderdir-m-2,j,k);

	if(cs==2)
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i,j-orderdir+m+2,k);

	if(cs==3)
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i,j+orderdir-m-2,k);

	if(cs==4)
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i-orderdir+m+2,j,k);

	if(cs==5)
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i,j,k+orderdir-m-2);

	if(cs==6)
	for(m=0;m<=orderdir-2;m++)
	y[m]=f(i,j,k-orderdir+m+2);

	y[orderdir-1]=wallvalue;

	if(ys==1 && dist<gamma*dx)
    {
    imagepoint(p,f,x_ip,val_ip,dist,cs);

    pos[orderdir-2] = x_ip;
    y[orderdir-2] = val_ip;
    }
		

	for(q=0; q<margin; ++q)
	{
	    y[orderdir+q]=0.0;

		for(m=0;m<orderdir;m++)
		{
			weight=1.0;
			for(n=0;n<orderdir;++n)
			{
			if(m!=n)
			weight*=(x[q]-pos[n])/(pos[m]-pos[n]+1.0e-20);
			}
		y[orderdir+q]+=weight*y[m];
		}
	}

// write extrapolated ghost cell values into f()

	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k) = y[orderdir+q-1+ys];

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k) = y[orderdir+q-1+ys];

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k) = y[orderdir+q-1+ys];

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k) = y[orderdir+q-1+ys];

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1) = y[orderdir+q-1+ys];

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1) = y[orderdir+q-1+ys];
}
