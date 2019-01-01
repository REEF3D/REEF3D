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

#include"pressure_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"

pressure_void::pressure_void(lexer* p)
{

}


pressure_void::~pressure_void()
{
}

void pressure_void::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    //LOOP
    //a->press()-
}

void pressure_void::ucorr(fdm* a,field& b,lexer*p)
{
}

void pressure_void::vcorr(fdm* a,field& b,lexer*p)
{
}

void pressure_void::wcorr(fdm* a,field& b,lexer*p)
{
}

void pressure_void::upgrad(lexer*p,fdm* a)
{
    if(p->D38==1)
    ULOOP
    {
	a->F(i,j,k)-=PORVAL1*fabs(p->W22)*(a->eta(i+1,j)-a->eta(i,j))/p->DXP[IP];
    }
    
    if(p->D38==2)
    ULOOP
    {
	a->F(i,j,k)-=PORVAL1*(a->eta(i+1,j)-a->eta(i,j))/p->DXP[IP];
    }
}

void pressure_void::vpgrad(lexer*p,fdm* a)
{
    if(p->D38==1)
    VLOOP
    {
	a->G(i,j,k)-=PORVAL2*fabs(p->W22)*(a->eta(i,j+1)-a->eta(i,j))/p->DYP[JP];
    }
    
    if(p->D38==2)
    VLOOP
    {
	a->G(i,j,k)-=PORVAL2*(a->eta(i,j+1)-a->eta(i,j))/p->DYP[JP];
    }
}

void pressure_void::wpgrad(lexer*p,fdm* a)
{
    double z1,z2;
    
    if(p->D38==1)
    WLOOP
    {
    z1 = p->ZP[KP];
    z2 = p->ZP[KP1];
	a->H(i,j,k)-=PORVAL3*fabs(p->W22)*(-z2+z1)/p->DZP[KP];
    }
    
    if(p->D38==2)
    WLOOP
    {
    z1 = p->ZP[KP];
    z2 = p->ZP[KP1];
	a->H(i,j,k)-=PORVAL3*(-z2+z1)/p->DZP[KP];
    }
}

void pressure_void::fillapu(lexer*p,fdm* a)
{
}

void pressure_void::fillapv(lexer*p,fdm* a)
{
}

void pressure_void::fillapw(lexer*p,fdm* a)
{
}

void pressure_void::rhs(lexer *p, fdm* a, ghostcell *pgc, field& uu, field& vv, field& ww, double alpha)
{
}

void pressure_void::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}


