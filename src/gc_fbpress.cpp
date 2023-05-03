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

void ghostcell::fbpress(lexer *p,field& f, double dist, int gcv, int bc, int cs)
{
	if(p->X33==1)
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
	f(i,j,k-q-1)=f(i,j,k);

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=f(i,j,k);
	}
    
		/*																		// Extrapolation such that dp/dx=0 at interface
		for(q=0;q<margin;++q)
		{
			double x1 = p->pos_z() - p->DXM;
			double x2 = p->pos_z();
			double xGamma = a->fb(i,j,k) + p->pos_z();
			double f1 = f(i,j,k-1);
			double f2 = f(i,j,k);
			double x3 = p->pos_z() + (q+1)*p->DXM;
			
			double a = (f2 - f1)/(x2*x2 - x1*x1 + 2.0*xGamma*(x1 - x2));
			
			f(i,j,k+q+1) = f1 + a*(x3*x3 - 2.0*xGamma*x3 - x1*x1 + 2.0*xGamma*x1);
		}*/
	
	if(p->X33==2)
	{
	double ui,vi,wi,un,vn,wn;
	
	ui = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
	vi = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
	wi = p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi;
	
	un = p->ufbn + (p->pos_z()-p->zgn)*p->qfbn - (p->pos_y()-p->ygn)*p->rfbn;
	vn = p->vfbn + (p->pos_x()-p->xgn)*p->rfbn - (p->pos_z()-p->zgn)*p->pfbn;
	wn = p->wfbn + (p->pos_y()-p->ygn)*p->pfbn - (p->pos_x()-p->xgn)*p->qfbn;

	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=f(i,j,k) - (p->DXM*(double(q)+1.0))*(ui-un)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=f(i,j,k) + (p->DXM*(double(q)+1.0))*(vi-vn)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=f(i,j,k) - (p->DXM*(double(q)+1.0))*(vi-vn)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=f(i,j,k) + (p->DXM*(double(q)+1.0))*(ui-un)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=f(i,j,k) - (p->DXM*(double(q)+1.0))*(wi-wn)/(p->dt>1.0e-20?p->dt:1.0e20);

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=f(i,j,k) + (p->DXM*(double(q)+1.0))*(wi-wn)/(p->dt>1.0e-20?p->dt:1.0e20);
	}
     
}
