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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"field.h"
#include"vec.h"
#include"fdm.h"

void ghostcell::gc_periodic(lexer *p, field& f, int gcv, int cs)
{
    double val1,val2,val3;
    
    if(cs==1)
    JLOOP
    KLOOP
    PCHECK
    {
    // 4 to 1 coupling
    i=p->knox-1;
    
    val1 = f(i,j,k);
    val2 = f(i-1,j,k);
    val3 = f(i-2,j,k);
    

    i=0;
    
	f(i-1,j,k) = val1;
    f(i-2,j,k) = val2;
    f(i-3,j,k) = val3;
    
    
    // 1 to 4 coupling
    i=0;
    
    val1 = f(i,j,k);
    val2 = f(i+1,j,k);
    val3 = f(i+2,j,k);
    
    i=p->knox-1;
    
	f(i+1,j,k) = val1;
    f(i+2,j,k) = val2;
    f(i+3,j,k) = val3;
    }
    
    
    if(cs==2)
    ILOOP
    KLOOP
    PCHECK
    {
    // 2 to 3 coupling
    j=p->knoy-1;
    
    val1 = f(i,j,k);
    val2 = f(i,j-1,k);
    val3 = f(i,j-2,k);
    
    j=0;
    
	f(i,j-1,k) = val1;
    f(i,j-2,k) = val2;
    f(i,j-3,k) = val3;
    
    // 3 to 2 coupling
    j=0;
    
    val1 = f(i,j,k);
    val2 = f(i,j+1,k);
    val3 = f(i,j+2,k);
    
    j=p->knoy-1;
    
	f(i,j+1,k) = val1;
    f(i,j+2,k) = val2;
    f(i,j+3,k) = val3;
    }
    
    
    if(cs==3)
    ILOOP
    JLOOP
    PCHECK
    {
    // 6 to 5 coupling
    k=p->knoz-1;
    
    val1 = f(i,j,k);
    val2 = f(i,j,k-1);
    val3 = f(i,j,k-2);
    
    k=0;
    
	f(i,j,k-1) = val1;
    f(i,j,k-2) = val2;
    f(i,j,k-3) = val3;
    
    // 5 to 6 coupling
    k=0;
    
    val1 = f(i,j,k);
    val2 = f(i,j,k+1);
    val3 = f(i,j,k+2);
    
    k=p->knoz-1;

	f(i,j,k+1) = val1;
    f(i,j,k+2) = val2;
    f(i,j,k+3) = val3;
    }
}

