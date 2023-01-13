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

#include"isusprhs.h"
#include"lexer.h"
#include"fdm.h"

isusprhs::isusprhs(lexer *p)
{
    d50=p->S20;
    ks=3.0*d50;
    gi=9.81;
    rhosed=2650.0;
    rhowat=p->W1;
    ws=1.1*(rhosed/rhowat-1.0)*gi*d50*d50;
}

isusprhs::~isusprhs()
{
}

void isusprhs::isuspsource(lexer* p,fdm* a,field& conc)
{
    count=0;
    LOOP
    {
	if(a->phi(i,j,k)>0.0)
	a->rhsvec.V[count]  += -ws*(conc(i,j,k+1)-conc(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
	
	++count;
    }
}

void isusprhs::sedfsf(lexer* p,fdm* a,field& conc)
{
    LOOP
    if(a->phi(i,j,k)<0.0)
    conc(i,j,k)=0.0;
}

void isusprhs::clearrhs(lexer* p, fdm* a)
{
    count=0;
    LOOP
    {
    a->rhsvec.V[count]=0.0;
	a->L(i,j,k)=0.0;
	++count;
    }
}
