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

#include"slice4a.h"
#include"lexer.h"

void slice4a::ggcpol(lexer* p)
{
    double val=0.0;
    int a,n,q,count;

	 for(n=0;n<p->ggcslcount4a;n++)
	 {
		 i=p->ggcsl4a[n][0];
		 j=p->ggcsl4a[n][1];

   
        val=0.0;
        count=0;
        
        iter=IJ;

        // 1
        for(q=1;q<=3;++q)
        if(p->gcslorig4a[p->mgcsl4a[iter]-10][0][q]==1)
        {
        val+=gcfeld[p->mgcsl4a[iter]-10][0][q];
        ++count;
        }
        
        // 4
        for(q=1;q<=3;++q)
        if(p->gcslorig4a[p->mgcsl4a[iter]-10][3][q]==1)
        {
        val+=gcfeld[p->mgcsl4a[iter]-10][3][q];
        ++count;
        }
        
        // 3
        for(q=1;q<=3;++q)
        if(p->gcslorig4a[p->mgcsl4a[iter]-10][2][q]==1)
        {
        val+=gcfeld[p->mgcsl4a[iter]-10][2][q];
        ++count;
        }
        
        // 2
        for(q=1;q<=3;++q)
        if(p->gcslorig4a[p->mgcsl4a[iter]-10][1][q]==1)
        {
        val+=gcfeld[p->mgcsl4a[iter]-10][1][q];
        ++count;
        }
		 
        if(count>0)
        val/=double(count);
    

        V[iter]=val; 
        
	 }
}









