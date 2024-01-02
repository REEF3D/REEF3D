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
Author: Tobias Martin
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

void sixdof_sflow::isource(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    {
    dfdx_plus = (press(i+1,j)-press(i,j))/p->DXP[IP];
    dfdx_min  = (press(i,j)-press(i-1,j))/p->DXP[IM1];
    
    dfdx = limiter(dfdx_plus,dfdx_min);
    
    dfdx = (press(i+1,j)-press(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
 
    d->F[IJK] += 1.0/p->W1*dfdx;
    }
    
    //SLICELOOP4
    //d->test2D(i,j) = press(i,j);
}

void sixdof_sflow::jsource(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    LOOP
    {
    dfdy_plus = (press(i,j+1)-press(i,j))/p->DYP[JP];
    dfdy_min  = (press(i,j)-press(i,j-1))/p->DYP[JM1];
    
    dfdy = limiter(dfdy_plus,dfdy_min);
    
    dfdy = (press(i,j+1)-press(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
        
    d->G[IJK] += 1.0/p->W1*dfdy;
    
    if(k==0)
    d->test2D(i,j) = 1.0/p->W1*dfdy;
    }
}

void sixdof_sflow::ksource(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}

void sixdof_sflow::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
	SLICELOOP1
    {
        b->F(i,j) += 1.0/p->W1*(press(i+1,j) - press(i,j))/p->DXP[IP];
    }
}

void sixdof_sflow::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
	SLICELOOP2
    {
        b->G(i,j) += 1.0/p->W1*(press(i,j+1) - press(i,j))/p->DYP[JP];
    }
    
    SLICELOOP4
    {
        b->test(i,j) = 1.0/p->W1*(press(i,j+1) - press(i,j))/p->DYP[JP];
    }
    pgc->gcsl_start4(p,b->test,50);
}

double sixdof_sflow::limiter(double v1, double v2)
{
    r=v2/(fabs(v1)>1.0e-10?v1:1.0e20);
    
    phival = (r*r + r)/(r*r+1.0);
    
    val = 0.5*phival*(v1+v2);

    return val;	
}
