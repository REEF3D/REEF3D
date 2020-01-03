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

#include"sflow_boussinesq_abbott.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"slice1.h"
#include"slice2.h"
 
sflow_boussinesq_abbott::sflow_boussinesq_abbott(lexer *p, fdm2D *b) : Pxx(p),Pxx_n(p),Qyy(p),Qyy_n(p)
{
}

sflow_boussinesq_abbott::~sflow_boussinesq_abbott()
{
}

void sflow_boussinesq_abbott::ini(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    psi1_calc(p,b,pgc,P,Q,b->eta,0.0);
    psi2_calc(p,b,pgc,P,Q,b->eta,0.0);
}

void sflow_boussinesq_abbott::psi1(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    if(p->count==0)
    ini(p,b,pgc,P,Q);
    
    SLICELOOP1
    Pxx_n(i,j) = Pxx(i,j);

    psi1_calc(p,b,pgc,P,Q,eta,alpha);
    
    // add to source
    if(p->count>1)
    SLICELOOP1
    {
    b->F(i,j) +=   ((1.0/3.0)*pow(b->depth(i,j),2.0)*(Pxx(i,j) - Pxx_n(i,j)))/(alpha*p->dt);
    }
}

void sflow_boussinesq_abbott::psi2(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    SLICELOOP2
    Qyy_n(i,j) = Qyy(i,j);

    psi2_calc(p,b,pgc,P,Q,eta,alpha);
    
    // add to source
    SLICELOOP2
    b->G(i,j) += (1.0/3.0)*pow(b->depth(i,j),2.0)*(Qyy(i,j) - Qyy_n(i,j))/(alpha*p->dt);
}

void sflow_boussinesq_abbott::psi1_calc(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    SLICELOOP1
    Pxx(i,j) = (P(i+1,j) - 2.0*P(i,j) + P(i-1,j))/(p->dx*p->dx);

    /*
    if(p->mpirank==0)
    {
    Pxx(0,0)  = 0.0;
    }*/
}

void sflow_boussinesq_abbott::psi2_calc(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    SLICELOOP2 
    Qyy(i,j) = (Q(i,j+1) - 2.0*Q(i,j) + Q(i,j-1))/(p->dx*p->dx);
}