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
#include"ghostcell.h"
#include"slice.h"

void iowave::eta_relax(lexer *p, ghostcell *pgc, slice &f)
{
	count=0;
    SLICELOOP4
    {
		dg = distgen(p);
		db = distbeach(p);
        

		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            f(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*etaval[count] + relax4_wg(i,j) * f(i,j);
            ++count;
            }
		}
        
		
		// Numerical Beach
		if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(p->A10!=3 || p->A348==1 || p->A348==2)
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            f(i,j) = relax4_nb(i,j)*f(i,j);
            }
        }
    }
}

void iowave::um_relax(lexer *p, ghostcell *pgc, slice &P, slice &bed, slice &eta)
{
 
    count=0;
    SLICELOOP1
    {
		dg = distgen(p);
		db = distbeach(p);
        
    
        // Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            P(i,j) = (1.0-relax1_wg(i,j))*ramp(p) * uval[count] + relax1_wg(i,j)*P(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            P(i,j) = relax1_nb(i,j)*P(i,j);
        }
    }
}

void iowave::vm_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
    count=0;
    SLICELOOP2
    {
		dg = distgen(p);
		db = distbeach(p);

        // Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            Q(i,j) = (1.0-relax1_wg(i,j))*ramp(p) * vval[count] + relax1_wg(i,j)*Q(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            Q(i,j) = relax1_nb(i,j)*Q(i,j);
        }
    }
}

void iowave::wm_relax(lexer *p, ghostcell *pgc, slice &W, slice &bed, slice &eta)
{
    
    count=0;
    SLICELOOP4
    {

		dg = distgen(p);
		db = distbeach(p);

        // Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            W(i,j) = (1.0-relax4_wg(i,j))*ramp(p) * wval[count] + relax4_wg(i,j)*W(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            W(i,j) = relax4_nb(i,j)*W(i,j);
        }
    }
}

void iowave::ws_relax(lexer *p, ghostcell *pgc, slice &W, slice &bed, slice &eta)
{
	double wval=0.0;
    
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

        z=eta(i,j);

        wval = wave_w(p,pgc,xg,yg,z);
        
        
        // Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            W(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*wval + relax4_wg(i,j)*W(i,j);
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            W(i,j) = relax4_nb(i,j)*W(i,j);
        }
    }
}

void iowave::pm_relax(lexer *p, ghostcell *pgc, slice &f)
{
	
    SLICELOOP4
    {
		xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2)
        {
            // Zone 1
            if(dg<1.0e20)
            f(i,j) = relax4_wg(i,j) * f(i,j);
		}
		
		// Numerical Beach
		if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            f(i,j) = relax4_nb(i,j)*f(i,j);
        }
    }
}
