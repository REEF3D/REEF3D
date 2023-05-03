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
#include"sflow_ediff.h"
#include"lexer.h"
#include"fdm2D.h"

sflow_ediff::sflow_ediff(lexer* p)
{
}

sflow_ediff::~sflow_ediff()
{
}

void sflow_ediff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{
    double visc=p->W2;
    
	SLICELOOP1
    {
	b->F(i,j) +=  ((visc+0.5*(b->eddyv(i,j) + b->eddyv(i+1,j)))/(p->DXM*p->DXM))*
    
                (2.0*(u(i+1,j) - 2.0*u(i,j) + u(i-1,j))
                    +(u(i,j+1) - 2.0*u(i,j) + u(i,j-1))
                
                + (v(i+1,j)-v(i,j)) - (v(i+1,j-1)-v(i,j-1)));
                                        
    }

}

void sflow_ediff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, double alpha)
{
    double visc=p->W2;
    
	SLICELOOP2
    {
	b->G(i,j) +=  ((visc+0.5*(b->eddyv(i,j) + b->eddyv(i,j+1)))/(p->DXM*p->DXM))*
    
                (     (v(i+1,j) - 2.0*v(i,j) + v(i-1,j))
                + 2.0*(v(i,j+1) - 2.0*v(i,j) + v(i,j-1))
                
                + (u(i,j+1)-u(i,j)) - (u(i-1,j+1)-u(i-1,j)));
    }
}

void sflow_ediff::diff_w(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &u, slice &v, slice &w, double alpha)
{
    double visc=p->W2;
    
	SLICELOOP4
    {
	b->L(i,j) +=  ((visc+0.5*(b->eddyv(i,j) + b->eddyv(i,j+1)))/(p->DXM*p->DXM))*
    
                (     (v(i+1,j) - 2.0*v(i,j) + v(i-1,j))
                + 2.0*(v(i,j+1) - 2.0*v(i,j) + v(i,j-1))
                
                + (u(i,j+1)-u(i,j)) - (u(i-1,j+1)-u(i-1,j)));
    }
}

void sflow_ediff::diff_scalar(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double sig, double alpha)
{
    
}
