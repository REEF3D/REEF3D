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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void iowave::U_relax(lexer *p, ghostcell *pgc, double *U)
{
    count=0;
    LOOP
    {
         dg = distgen(p);
		db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            U[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*uval[count] + relax4_wg(i,j)*U[IJK];
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            U[IJK] = relax4_nb(i,j)*U[IJK];
        }
    }
}

void iowave::V_relax(lexer *p, ghostcell *pgc, double *V)
{
    count=0;
    if(p->j_dir==0)
    LOOP
    {
        dg = distgen(p);
		db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            V[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*vval[count] + relax4_wg(i,j)*V[IJK];
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1||p->B99==2||beach_relax==1)
		{	
            // Zone 2
            if(db<1.0e20)
            V[IJK] = relax4_nb(i,j)*V[IJK];
        }
    }
}

void iowave::W_relax(lexer *p, ghostcell *pgc, double *W)
{
    count=0;
    LOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            W[IJK] = (1.0-relax4_wg(i,j)) * ramp(p)* wval[count] + relax4_wg(i,j)*W[IJK];
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            W[IJK] = relax4_nb(i,j)*W[IJK];
        }
    }		
}

void iowave::P_relax(lexer *p, ghostcell *pgc, double *P)
{
    FLOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
        // Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
        {            
            // Zone 2
            if(db<1.0e20)

            P[FIJK] = relax4_nb(i,j)*P[FIJK];
        }
    }	
}


