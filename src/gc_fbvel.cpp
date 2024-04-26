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

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"

void ghostcell::fbvel1(lexer *p,field& f, double dist, int gcv, int bc, int cs)
{
	double vel,velpar;
	
	vel = p->ufbi + (p->pos_z()-p->zg)*p->qfbi - (p->pos_y()-p->yg)*p->rfbi;
	
	if(p->X31==1)
	velpar=f(i,j,k);
	
	if(p->X31==2)
	velpar=0.0;
	
	if(p->X31==3)
	velpar=vel;
	
	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=velpar;

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=velpar;

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=velpar;

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=velpar;

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=velpar;

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=velpar;
	
	if(p->X31==4)
	{
		for(q = 0; q<margin; ++q)
		{			
			if(cs==1)	// right side
			{	
				if ((a->fb(i,j,k) - p->DXN[IP]/2.0) > 0.0)
				{
					double delta_x = (a->fb(i,j,k) - p->DXN[IP]/2.0)/(a->fb(i,j,k) + p->DXN[IP]/2.0);
					
					f(i-q-1, j, k) = f(i,j,k)*delta_x + vel*(1.0 - delta_x);
				}
				else
				{
					f(i-q-1, j, k) = vel;
				}
			}
				
			if(cs==2)	// front side
			{				
				// u(i,j+q+1,k) is always in solid
				f(i, j+q+1, k) = vel;
				
				// u(i,j,k)  always last velocity in fluid and at same height as pressure
				double delta_y = a->fb(i,j,k)/a->fb(i,j-1,k);
				
				f(i, j, k) = f(i,j-1,k)*delta_y + vel*(1.0 - delta_y);
			}
				 
			if(cs==3)	// back side
			{
				// u(i,j-q-1,k) is always in solid
				f(i, j-q-1, k) = vel;
				
				// u(i,j,k)  always last velocity in fluid and at same height as pressure
				double delta_y = a->fb(i,j,k)/a->fb(i,j+1,k);
				
				f(i, j, k) = f(i,j+1,k)*delta_y + vel*(1.0 - delta_y);
			}

			if(cs==4) // left side
			{			
				int iq1 = IP1;
				if (q==1) iq1 = IP2;
				if (q==2) iq1 = IP3;			
						
				if ((a->fb(i+q+1,j,k) - p->DXN[iq1]/2.0) > 0.0)
				{
					double delta_x = (a->fb(i+q+1,j,k) - p->DXN[iq1]/2.0)/(a->fb(i,j,k) - p->DXN[IP]/2.0);
					
					f(i+q+1, j, k) = f(i,j,k)*delta_x + vel*(1.0 - delta_x);
				}
				else
				{
					f(i+q+1, j, k) = vel;
				}
			}
				
			if(cs==5) // top 
			{			
				// u(i,j,k-q-1) is always in solid
				f(i, j, k-q-1) = vel;
				
				// u(i,j,k)  always last velocity in fluid and at same height as pressure
				double delta_z = a->fb(i,j,k)/a->fb(i,j,k+1);
				
				f(i, j, k) = f(i,j,k+1)*delta_z + vel*(1.0 - delta_z);
			}
			if(cs==6)	// bottom 
			{	
				// u(i,j,k+q+1) is always in solid
				f(i, j, k+q+1) = vel;
				
				// u(i,j,k)  always last velocity in fluid and at same height as pressure
				double delta_z = a->fb(i,j,k)/a->fb(i,j,k-1);
				
				f(i, j, k) = f(i,j,k-1)*delta_z + vel*(1.0 - delta_z);
			}
		}
	}


	if(p->X34==1)
    {
        if(cs==1&&bc==43)
        f(i-1,j,k)=0.5*(vel+f(i,j,k));

        if(cs==2&&bc==43)
        f(i,j+1,k)=0.5*(velpar+f(i,j,k));

        if(cs==3&&bc==43)
        f(i,j-1,k)=0.5*(velpar+f(i,j,k));

        if(cs==4&&bc==43)
        f(i+1,j,k)=0.5*(vel+f(i,j,k));

        if(cs==5&&bc==43)
        f(i,j,k-1)=0.5*(velpar+f(i,j,k));

        if(cs==6&&bc==43)
        f(i,j,k+1)=0.5*(velpar+f(i,j,k));
    }
}

void ghostcell::fbvel2(lexer *p,field& f, double dist, int gcv, int bc, int cs)
{
	double vel,velpar;
	
	vel = p->vfbi + (p->pos_x()-p->xg)*p->rfbi - (p->pos_z()-p->zg)*p->pfbi;
	
	if(p->X31==1)
	velpar=f(i,j,k);
	
	if(p->X31==2)
	velpar=0.0;
	
	if(p->X31==3)
	velpar=vel;
	
	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=velpar;
	
	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=velpar;

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=velpar;

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=velpar;

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=velpar;

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=velpar;
    
	
	if(p->X31==4)
	{
		for(q = 0; q<margin; ++q)
		{			
			if(cs==1)	// right side
			{
				// v(i-q-1,j,k) is always in solid
				f(i-q-1, j, k) = vel;
				
				// v(i,j,k) always last velocity in fluid and at same height as pressure
				double delta_x = a->fb(i,j,k)/a->fb(i+1,j,k);
				
				f(i, j, k) = f(i+1,j,k)*delta_x + vel*(1.0 - delta_x);
				
			}
				
			if(cs==2)	// front side --> cs=4 of u
			{				
				int jq1 = JP1;
				if (q==1) jq1 = JP2;
				if (q==2) jq1 = JP3;			
						
				if ((a->fb(i,j+q+1,k) - p->DYN[jq1]/2.0) > 0.0)
				{
					double delta_y = (a->fb(i,j+q+1,k) - p->DYN[jq1]/2.0)/(a->fb(i,j,k) - p->DYN[JP]/2.0);
					
					f(i, j+q+1, k) = f(i,j,k)*delta_y + vel*(1.0 - delta_y);
				}
				else
				{
					f(i, j+q+1, k) = vel;
				}
			}
				
			if(cs==3)	// back side --> cs=1 of u
			{
				if ((a->fb(i,j,k) - p->DYN[JP]/2.0) > 0.0)
				{
					double delta_y = (a->fb(i,j,k) - p->DYN[JP]/2.0)/(a->fb(i,j,k) + p->DYN[JP]/2.0);
					
					f(i, j-q-1, k) = f(i,j,k)*delta_y + vel*(1.0 - delta_y);
				}
				else
				{
					f(i, j-q-1, k) = vel;
				}
			}

			if(cs==4) // left side
			{
				// v(i+q+1,j,k) is always in solid
				f(i+q+1, j, k) = vel;
				
				// v(i,j,k) always last velocity in fluid and at same height as pressure
				double delta_x = a->fb(i,j,k)/a->fb(i-1,j,k);
				
				f(i, j, k) = f(i-1,j,k)*delta_x + vel*(1.0 - delta_x);
			}
				
			if(cs==5) // top 
			{			
				// v(i,j,k-q-1) is always in solid
				f(i, j, k-q-1) = vel;
				
				// v(i,j,k)  always last velocity in fluid and at same height as pressure
				double delta_z = a->fb(i,j,k)/a->fb(i,j,k+1);
				
				f(i, j, k) = f(i,j,k+1)*delta_z + vel*(1.0 - delta_z);
			}
			if(cs==6)	// bottom 
			{	
				// v(i,j,k+q+1) is always in solid
				f(i, j, k+q+1) = vel;
				
				// v(i,j,k)  always last velocity in fluid and at same height as pressure
				double delta_z = a->fb(i,j,k)/a->fb(i,j,k-1);
				
				f(i, j, k) = f(i,j,k-1)*delta_z + vel*(1.0 - delta_z);
			}
		}
	}	
	
	
    if(p->X34==1)
    {
        if(cs==1&&bc==43)
        f(i-1,j,k)=0.5*(velpar+f(i,j,k));

        if(cs==2&&bc==43)
        f(i,j+1,k)=0.5*(vel+f(i,j,k));

        if(cs==3&&bc==43)
        f(i,j-1,k)=0.5*(vel+f(i,j,k));

        if(cs==4&&bc==43)
        f(i+1,j,k)=0.5*(vel+f(i,j,k));

        if(cs==5&&bc==43)
        f(i,j,k-1)=0.5*(velpar+f(i,j,k));

        if(cs==6&&bc==43)
        f(i,j,k+1)=0.5*(velpar+f(i,j,k));
    }
}

void ghostcell::fbvel3(lexer *p,field& f, double dist, int gcv, int bc, int cs)
{	
	double vel,velpar;
	
	vel = p->wfbi + (p->pos_y()-p->yg)*p->pfbi - (p->pos_x()-p->xg)*p->qfbi;
	
	if(p->X31==1)
	velpar=f(i,j,k);
	
	if(p->X31==2)
	velpar=0.0;
	
	if(p->X31==3)
	velpar=vel;

	if(cs==1)
	for(q=0;q<margin;++q)
	f(i-q-1,j,k)=velpar;

	if(cs==2)
	for(q=0;q<margin;++q)
	f(i,j+q+1,k)=velpar;

	if(cs==3)
	for(q=0;q<margin;++q)
	f(i,j-q-1,k)=velpar;

	if(cs==4)
	for(q=0;q<margin;++q)
	f(i+q+1,j,k)=velpar;

	if(cs==5)
	for(q=0;q<margin;++q)
	f(i,j,k-q-1)=velpar;

	if(cs==6)
	for(q=0;q<margin;++q)
	f(i,j,k+q+1)=velpar;
	
	
	if(p->X31==4)
	{
		for(q = 0; q<margin; ++q)
		{		
			if(cs==1)	// right side										
			{
				// w(i-q-1,j,k) is always in solid
				f(i-q-1, j, k) = vel;
				
				// w(i,j,k) always last velocity in fluid and at same height as pressure
				double delta_x = a->fb(i,j,k)/a->fb(i+1,j,k);
				
				f(i, j, k) = f(i+1,j,k)*delta_x + vel*(1.0 - delta_x);
			}
				
			if(cs==2)	// front side								
			{	
				// w(i,j+q+1,k) is always in solid
				f(i, j+q+1, k) = vel;
				
				// w(i,j,k) always last velocity in fluid and at same height as pressure
				double delta_y = a->fb(i,j,k)/a->fb(i,j-1,k);
				
				f(i, j, k) = f(i,j-1,k)*delta_y + vel*(1.0 - delta_y);
			}
				
			if(cs==3)	// back side											
			{
				// w(i-q-1,j,k) is always in solid
				f(i, j-q-1, k) = vel;
				
				// w(i,j,k) always last velocity in fluid and at same height as pressure
				double delta_y = a->fb(i,j,k)/a->fb(i,j+1,k);
				
				f(i, j, k) = f(i,j+1,k)*delta_y + vel*(1.0 - delta_y);
			}

			if(cs==4) // left side											
			{
				// w(i+q+1,j,k) is always in solid
				f(i+q+1, j, k) = vel;
				
				// w(i,j,k) always last velocity in fluid and at same height as pressure
				double delta_x = a->fb(i,j,k)/a->fb(i-1,j,k);
				
				f(i, j, k) = f(i-1,j,k)*delta_x + vel*(1.0 - delta_x);
			}	
				
			if(cs==5) // top --> cs=1 for u
			{	
				if ((a->fb(i,j,k) - p->DZN[KP]/2.0) > 0.0)
				{
					double delta_z = (a->fb(i,j,k) - p->DZN[KP]/2.0)/(a->fb(i,j,k) + p->DZN[KP]/2.0);
					
					f(i, j, k-q-1) = f(i,j,k)*delta_z + vel*(1.0 - delta_z);
				}
				else
				{
					f(i, j, k-q-1) = vel;
				}			

			}
			if(cs==6)	// bottom --> cs=4 for u
			{	
				int kq1 = KP1;
				if (q==1) kq1 = KP2;
				if (q==2) kq1 = KP3;			
						
				if ((a->fb(i,j,k+q+1) - p->DZN[kq1]/2.0) > 0.0)
				{
					double delta_z = (a->fb(i,j,k+q+1) - p->DZN[kq1]/2.0)/(a->fb(i,j,k) - p->DZN[KP]/2.0);
					
					f(i, j, k+q+1) = f(i,j,k)*delta_z + vel*(1.0 - delta_z);
				}
				else
				{
					f(i, j, k+q+1) = vel;
				}				
			}
		}
	}	
	
	
	if(p->X34==1)
    {
        if(cs==1&&bc==43)
        f(i-1,j,k)=0.5*(velpar+f(i,j,k));

        if(cs==2&&bc==43)
        f(i,j+1,k)=0.5*(velpar+f(i,j,k));

        if(cs==3&&bc==43)
        f(i,j-1,k)=0.5*(velpar+f(i,j,k));

        if(cs==4&&bc==43)
        f(i+1,j,k)=0.5*(velpar+f(i,j,k));

        if(cs==5&&bc==43)
        f(i,j,k-1)=0.5*(vel+f(i,j,k));

        if(cs==6&&bc==43)
        f(i,j,k+1)=0.5*(vel+f(i,j,k));
	}
	
}


