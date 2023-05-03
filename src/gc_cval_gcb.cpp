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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"

void ghostcell::cval_gcb4(lexer* p, fdm* a, fieldint &cval)
{
	GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	
	p->gcb4[n][5]=cval(i,j,k);
	}
    
    
    if(p->periodic1==1)
    {
        i=0;
        n=0;
        JLOOP
        KLOOP
        PCHECK
        {
        p->gc4periodic[0][n]=cval(i,j,k);
        ++n;
        }
        p->gc4periodic_count[0]=n;
        
        i=p->knox-1;
        n=0;
        JLOOP
        KLOOP
        PCHECK
        {
        p->gc4periodic[3][n]=cval(i,j,k);
        ++n;
        }
        p->gc4periodic_count[3]=n;
    }
    
    if(p->periodic2==1)
    {
        j=0;
        n=0;
        ILOOP
        KLOOP
        PCHECK
        {
        p->gc4periodic[2][n]=cval(i,j,k);
        ++n;
        }
        p->gc4periodic_count[2]=n;
        
        j=p->knoy-1;
        n=0;
        ILOOP
        KLOOP
        PCHECK
        {
        p->gc4periodic[1][n]=cval(i,j,k);
        ++n;
        }
        p->gc4periodic_count[1]=n;
    }
    
    if(p->periodic3==1)
    {
        k=0;
        n=0;
        ILOOP
        JLOOP
        PCHECK
        {
        p->gc4periodic[4][n]=cval(i,j,k);
        ++n;
        }
        p->gc4periodic_count[4]=n;
        
        k=p->knoz-1;
        n=0;
        ILOOP
        JLOOP
        PCHECK
        {
        p->gc4periodic[5][n]=cval(i,j,k);
        ++n;
        }
        p->gc4periodic_count[5]=n;
    }
}

void ghostcell::cval_gcb4a(lexer* p, fdm* a, fieldint &cval)
{
	GC4ALOOP
    {
    i=p->gcb4a[n][0];
    j=p->gcb4a[n][1];
    k=p->gcb4a[n][2];
	
	p->gcb4a[n][5]=cval(i,j,k);
	}
    
    
    if(p->periodic1==1)
    {
        i=0;
        n=0;
        JLOOP
        KLOOP
        PSOLIDCHECK
        {
        p->gc4aperiodic[0][n]=cval(i,j,k);
        ++n;
        }
        p->gc4aperiodic_count[0]=n;
        
        i=p->knox-1;
        n=0;
        JLOOP
        KLOOP
        PSOLIDCHECK
        {
        p->gc4aperiodic[3][n]=cval(i,j,k);
        ++n;
        }
        p->gc4aperiodic_count[3]=n;
    }
    
    if(p->periodic2==1)
    {
        j=0;
        n=0;
        ILOOP
        KLOOP
        PSOLIDCHECK
        {
        p->gc4aperiodic[2][n]=cval(i,j,k);
        ++n;
        }
        p->gc4aperiodic_count[2]=n;
        
        j=p->knoy-1;
        n=0;
        ILOOP
        KLOOP
        PSOLIDCHECK
        {
        p->gc4aperiodic[1][n]=cval(i,j,k);
        ++n;
        }
        p->gc4aperiodic_count[1]=n;
    }
    
    if(p->periodic3==1)
    {
        k=0;
        n=0;
        ILOOP
        JLOOP
        PSOLIDCHECK
        {
        p->gc4aperiodic[4][n]=cval(i,j,k);
        ++n;
        }
        p->gc4aperiodic_count[4]=n;
        
        k=p->knoz-1;
        n=0;
        ILOOP
        JLOOP
        PSOLIDCHECK
        {
        p->gc4aperiodic[5][n]=cval(i,j,k);
        ++n;
        }
        p->gc4aperiodic_count[5]=n;
    }
}

void ghostcell::cval_gcb6(lexer* p, fdm* a, fieldint &cval)
{
    int count1,count2;
    count1=count2=0;
    
	GC6LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	
	p->gcb6[n]=cval(i,j,k);
    
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
        p->gcin6[count1][3] = p->gcb6[n];
        ++count1;
        }
        
        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
        p->gcout6[count2][3] = p->gcb6[n];
        ++count2;
        }
	}
}

