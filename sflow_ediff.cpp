/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"sflow_ediff.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_ediff::sflow_ediff(lexer* p)
{
}

sflow_ediff::~sflow_ediff()
{
}

void sflow_ediff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double alpha)
{
    double visc=p->W2;
    double dudx,dvdy,dudy,dvdx;
    
    
	SLICELOOP1
    {
    if(p->A260==1)
    {
    dudx = (b->P(i+1,j) - b->P(i-1,j))/(2.0*p->dx);
    dvdy = (0.5*(b->Q(i,j)+b->Q(i-1,j)) - 0.5*(b->Q(i,j-1)+b->Q(i-1,j-1)))/(p->dx);
    dudy = (b->P(i,j+1) - b->P(i,j-1))/(2.0*p->dx);
    dvdx = (0.5*(b->Q(i+1,j)+b->Q(i+1,j-1)) - 0.5*(b->Q(i,j)+b->Q(i,j-1)))/(p->dx);
    
    visc = 4.0*sqrt(2.0*pow(dudx,2.0) + 2.0*pow(dvdy,2.0) + pow(dudy+dvdx,2.0))+ p->W2;
    }
    
	b->F(i,j) +=  (visc/(p->dx*p->dx))*(f(i+1,j) - 2.0*f(i,j) + f(i-1,j)
									    + f(i,j+1) - 2.0*f(i,j) + f(i,j-1));
                                        
    }

}

void sflow_ediff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double alpha)
{
    double visc=p->W2;
    double dudx,dvdy,dudy,dvdx;
    
	SLICELOOP2
    {
        if(p->A260==1)
        {
        dudx = (0.5*(b->P(i+1,j+1)+b->P(i+1,j)) - 0.5*(b->P(i,j)+b->P(i,j)))/(p->dx);
        dvdy = (b->Q(i,j+1) - b->Q(i,j-1))/(2.0*p->dx);
        dudy = (0.5*(b->P(i,j+1)+b->P(i-1,j+1)) - 0.5*(b->P(i,j)+b->P(i-1,j)))/(p->dx);
        dvdx = (b->Q(i+1,j) - b->Q(i-1,j))/(2.0*p->dx);
        
        visc = 4.0*sqrt(2.0*pow(dudx,2.0) + 2.0*pow(dvdy,2.0) + pow(dudy+dvdx,2.0))+ p->W2;
        }
        
	b->G(i,j) +=  (visc/(p->dx*p->dx))*(f(i+1,j) - 2.0*f(i,j) + f(i-1,j)
									    + f(i,j+1) - 2.0*f(i,j) + f(i,j-1));
    }
}