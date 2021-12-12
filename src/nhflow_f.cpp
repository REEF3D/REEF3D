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

#include"nhflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

nhflow_f::nhflow_f(lexer *p, fdm *a, ghostcell *pgc) 
{
    margin=3;
}

nhflow_f::~nhflow_f()
{
}


void nhflow_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}

void nhflow_f::kinematic_fsf(lexer *p, fdm *a, field &u, field &v, field &w, slice &eta, slice &eta_n, double alpha)
{
    double wval;
    
    GC4LOOP
    if(p->gcb4[n][3]==6 && p->gcb4[n][4]==3)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
	wval = (a->eta(i,j) - a->eta_n(i,j))/(p->dt)
    
         + 0.5*(u(i,j,k)+u(i-1,j,k))*((eta(i+1,j)-eta(i-1,j))/(p->DXP[IP]+p->DXP[IP1]))
    
         + 0.5*(v(i,j,k)+v(i,j-1,k))*((eta(i,j+1)-eta(i,j-1))/(p->DYP[JP]+p->DYP[JP1]));
    
	for(q=0;q<margin;++q)
	w(i,j,k+q) = wval; 
    }
    

    GC4LOOP
    if(p->gcb4[n][3]==5 && p->gcb4[n][4]==21)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    wval = - 0.5*(u(i,j,k)+u(i-1,j,k))*((a->depth(i+1,j)-a->depth(i-1,j))/(p->DXP[IP]+p->DXP[IP1]))
    
           - 0.5*(v(i,j,k)+v(i,j-1,k))*((a->depth(i,j+1)-a->depth(i,j-1))/(p->DYP[JP]+p->DYP[JP1]));
    
   /* for(q=0;q<margin;++q)
    {
	a->test(i,j,k-q) = wval;
    }
    */
    
    //wval = 0.0;
    
	for(q=0;q<margin;++q)
	w(i,j,k-q-1) = wval;
    }
}


