/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"grid.h"
#include"lexer.h"
#include"fdm.h"
#include"fieldint6.h"

void grid::column_pt1_update(lexer* p, cpt &C1)
{
    fieldint1 cval1(p);
    cval_update1(p,cval1);	
	column_pt1_assign(p,cval1,C1);
    cval_gcb1(p,cval1);
    cval_gcpara1(p,cval1);
}

void grid::column_pt2_update(lexer* p, cpt &C2)
{
    fieldint2 cval2(p);
    cval_update2(p,cval2);	
	column_pt2_assign(p,cval2,C2);
    cval_gcb2(p,cval2);
    cval_gcpara2(p,cval2);
}

void grid::column_pt3_update(lexer* p, cpt &C3)
{
    fieldint3 cval3(p);
    cval_update3(p,cval3);	
	column_pt3_assign(p,cval3,C3);
    cval_gcb3(p,cval3);
    cval_gcpara3(p,cval3);
}

void grid::column_pt4_update(lexer* p, cpt &C4)
{
    fieldint4 cval4(p);
    cval_update4(p,cval4);	
	column_pt4_assign(p,cval4,C4);
    cval_gcb4(p,cval4);
    cval_gcpara4(p,cval4);
}

void grid::column_pt4a_update(lexer* p, cpt &C4a)
{
    fieldint4a cval4a(p);
    cval_update4a(p,cval4a);	
	column_pt4a_assign(p,cval4a,C4a);
    cval_gcb4a(p,cval4a);
    cval_gcpara4a(p,cval4a);
}

void grid::column_pt6_update(lexer* p, cpt &C6)
{
    fieldint6 cval6(p);
    cval_update6(p,cval6);	
	column_pt6_assign(p,cval6,C6);
    cval_gcb6(p,cval6);
    cval_gcpara6(p,cval6);
}
