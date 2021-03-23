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

#include"bedload_noneq.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

bedload_noneq::bedload_noneq(lexer *p) : q0(p)
{
    /*
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    visc=p->W2;
    kappa=0.4;
    ks=p->S21*d50;
    Rstar=(rhosed-rhowat)/rhowat;
    Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);*/
}

bedload_noneq::~bedload_noneq()
{
}

void bedload_noneq::ini(lexer* p, fdm* a, ghostcell* pgc, slice &q)
{

}

void bedload_noneq::start(lexer* p, fdm* a, ghostcell* pgc, slice &q)
{

}
