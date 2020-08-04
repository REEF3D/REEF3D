/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"sflow_hydrostatic.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver2D.h"
#include"momentum.h"
#include"ioflow.h"
 
sflow_hydrostatic::sflow_hydrostatic(lexer* p, fdm2D *b)
{
}

sflow_hydrostatic::~sflow_hydrostatic()
{
}

void sflow_hydrostatic::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, slice &P, slice &Q, slice &Pn, slice &Qn, slice &ws, slice &eta, double alpha)
{
    
}

void sflow_hydrostatic::ucorr(lexer* p, fdm2D* b, slice& uvel, slice &eta,double alpha)
{	
}

void sflow_hydrostatic::vcorr(lexer* p, fdm2D* b, slice& vvel, slice &eta,double alpha)
{	
}

void sflow_hydrostatic::wcorr(lexer* p, fdm2D* b,double alpha, slice &uvel, slice &vvel, slice &ws)
{	
}

void sflow_hydrostatic::wcalc(lexer* p, fdm2D* b,double alpha, slice &uvel, slice &vvel, slice &ws)
{	
}

void sflow_hydrostatic::upgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
	if(p->A221>=1)
    SLICELOOP1
    b->F(i,j) -= fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) 
                                 - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM); 
}

void sflow_hydrostatic::vpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
	if(p->A221>=1)
        SLICELOOP2
        b->G(i,j) -= fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) 
                                 - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM);
}

void sflow_hydrostatic::wpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{	    
}





