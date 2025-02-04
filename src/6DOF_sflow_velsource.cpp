/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_nhf.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sixdof_sflow::isource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::jsource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::ksource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(p->X10==3)
	SLICELOOP1
    {
    dfdx = (press(i+1,j)-press(i,j))/(p->DXP[IP]);
    
    b->F(i,j) += dfdx/p->W1;
    }
    
    SLICELOOP4
    b->test(i,j) = press(i,j);
}

void sixdof_sflow::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(p->X10==3)
	SLICELOOP2
    {
    dfdy = (press(i,j+1)-press(i,j))/(p->DYP[JP]);
    
    b->G(i,j) += dfdy/p->W1;
    }
}

void sixdof_sflow::isource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_sflow::jsource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

void sixdof_sflow::ksource(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
}

double sixdof_sflow::limiter(double v1, double v2)
{
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);
    
    phival = (r*r + r)/(r*r+1.0);
    
    val = 0.5*phival*(v1+v2);

    return val;	
}
