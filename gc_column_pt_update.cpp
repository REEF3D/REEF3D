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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fieldint6.h"

void ghostcell::column_pt1_update(lexer* p, fdm* a)
{
    fieldint1 cval1(p);
    cvalupdate1(p,a,cval1);	
	column_pt1(p,a,cval1);
    cval_gcb1(p,a,cval1);
    cval_gcpara1(p,a,cval1);
}

void ghostcell::column_pt2_update(lexer* p, fdm* a)
{
    fieldint2 cval2(p);
    cvalupdate2(p,a,cval2);	
	column_pt2(p,a,cval2);
    cval_gcb2(p,a,cval2);
    cval_gcpara2(p,a,cval2);
}

void ghostcell::column_pt3_update(lexer* p, fdm* a)
{
    fieldint3 cval3(p);
    cvalupdate3(p,a,cval3);	
	column_pt3(p,a,cval3);
    cval_gcb3(p,a,cval3);
    cval_gcpara3(p,a,cval3);
}

void ghostcell::column_pt4_update(lexer* p, fdm* a)
{
    fieldint4 cval4(p);
    cvalupdate4(p,a,cval4);	
	column_pt4(p,a,cval4);
    cval_gcb4(p,a,cval4);
    cval_gcpara4(p,a,cval4);
}

void ghostcell::column_pt4a_update(lexer* p, fdm* a)
{
    fieldint4a cval4a(p);
    cvalupdate4a(p,a,cval4a);	
	column_pt4a(p,a,cval4a);
    cval_gcb4a(p,a,cval4a);
    cval_gcpara4a(p,a,cval4a);
}

void ghostcell::column_pt6_update(lexer* p, fdm* a)
{
    fieldint6 cval6(p);
    cvalupdate6(p,a,cval6);	
	column_pt6(p,a,cval6);
    cval_gcb6(p,a,cval6);
    cval_gcpara6(p,a,cval6);
}