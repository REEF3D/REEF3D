/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

void iowave::eta_relax(lexer *p, ghostcell *pgc, slice &f)
{
	count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        
		
		// Wave Generation
        if(p->B98==1 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            f(i,j) = ra1(p,dg)*etaval[count] + (1.0-ra1(p,dg))*f(i,j);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            f(i,j) = ra2(p,dg)*etaval[count] + (1.0-ra2(p,dg))*f(i,j);
            ++count;
            }
		}

        if(p->B98==2 && h_switch==1)
        {
            
            // Zone 1
            if(dg<dist1)
            {
            //cout<<"DG: "<<dg<<" XG: "<<xg<<" rb1: "<<rb1(p,dg)<<endl;
            f(i,j) = (1.0-rb1(p,dg))*ramp(p)*etaval[count] + rb1(p,dg) * f(i,j);
            ++count;
            }
		}
        
		
		// Numerical Beach
		if(p->B99==1)
		{
		      // Zone 3
            if(db<dist3)
            f(i,j) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*f(i,j);
        }
		
		if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            f(i,j) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*f(i,j);
        }
    }
}

void iowave::um_relax(lexer *p, ghostcell *pgc, slice &P, slice &bed, slice &eta)
{
 
    count=0;
    SLICELOOP1
    {
        xg = xgen1(p);
        yg = ygen1(p);
		dg = distgen(p);
		db = distbeach(p);
        
    
        // Wave Generation
        if(p->B98==1 && u_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            P(i,j) = ra1(p,dg) * uval[count] + (1.0-ra1(p,dg))*P(i,j);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            P(i,j) = ra2(p,dg) * uval[count] + (1.0-ra2(p,dg))*P(i,j);
            ++count;
            }
		}
		
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            P(i,j) = (1.0-rb1(p,dg))*ramp(p) * uval[count] + rb1(p,dg)*P(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            P(i,j) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*P(i,j);
        }

        if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            P(i,j) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*P(i,j);
        }
    }
}

void iowave::vm_relax(lexer *p, ghostcell *pgc, slice &Q, slice &bed, slice &eta)
{
    count=0;
    SLICELOOP2
    {
        xg = xgen2(p);
        yg = ygen2(p);
		dg = distgen(p);
		db = distbeach(p);

        // Wave Generation
        if(p->B98==1 && v_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            Q(i,j) = ra1(p,dg) * vval[count]  + (1.0-ra1(p,dg))*Q(i,j);
            ++count;
            }

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            Q(i,j) = ra2(p,dg) * vval[count] + (1.0-ra2(p,dg))*Q(i,j);
            ++count;
            }
			
		}
		
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            Q(i,j) = (1.0-rb1(p,dg))*ramp(p) * vval[count] + rb1(p,dg)*Q(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            Q(i,j) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*Q(i,j);
        }

        if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            Q(i,j) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*Q(i,j);
        }
    }
}

void iowave::wm_relax(lexer *p, ghostcell *pgc, slice &W, slice &bed, slice &eta)
{
    
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

        // Wave Generation
        if(p->B98==1 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            W(i,j) = ra1(p,dg) * wval[count]  + (1.0-ra1(p,dg))*W(i,j);
            ++count;
            }
            

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            {
            W(i,j) = ra2(p,dg) * wval[count] + (1.0-ra2(p,dg))*W(i,j);
            ++count;
            }
			
		}
		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
            W(i,j) = (1.0-rb1(p,dg))*ramp(p) * wval[count] + rb1(p,dg)*W(i,j);
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            W(i,j) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*W(i,j);
        }

        if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            W(i,j) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*W(i,j);
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
        if(p->B98==1 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            W(i,j) = ra1(p,dg) * wval  + (1.0-ra1(p,dg))*W(i,j);

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            W(i,j) = ra2(p,dg) * wval + (1.0-ra2(p,dg))*W(i,j);
			
		}
		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            W(i,j) = (1.0-rb1(p,dg))*ramp(p)*wval + rb1(p,dg)*W(i,j);
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
            // Zone 3
            if(db<dist3)
            W(i,j) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*W(i,j);
        }

        if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            W(i,j) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*W(i,j);
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
        if(p->B98==1)
        {
            // Zone 1
            if(dg<dist1)
            f(i,j) = ra1(p,dg)*0.0 + (1.0-ra1(p,dg))*f(i,j);

            // Zone 2
            if(dg>=dist1 && dg<dist2)
            f(i,j) = ra2(p,dg)*0.0 + (1.0-ra2(p,dg))*f(i,j);
		}

        if(p->B98==2)
        {
            // Zone 1
            if(dg<dist1)
            f(i,j) = (1.0-rb1(p,dg))*0.0 + rb1(p,dg) * f(i,j);
		}
		
		// Numerical Beach
		if(p->B99==1)
		{
		      // Zone 3
            if(db<dist3)
            f(i,j) = (1.0-ra3(p,db))*0.0 + ra3(p,db)*f(i,j);
        }
		
		if(p->B99==2)
		{
            // Zone 3
            if(db<dist3)
            f(i,j) = (1.0-rb3(p,db))*0.0 + rb3(p,db)*f(i,j);
        }
    }
}