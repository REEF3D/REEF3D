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

#include"field1.h"
#include"cart1.h"
#include"lexer.h"

void field1::ggcpol(lexer* p)
{
    double val=0.0;
    int a,n,q,count;

	 for(n=0;n<cart1::ggccount;n++)
	 {
		 i=cart1::ggc[n][0];
		 j=cart1::ggc[n][1];
		 k=cart1::ggc[n][2];

		 iter=(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin;
		 val=0.0;
        count=0;

        // 1
        for(q=1;q<=3;++q)
        if(p->gcorig1[pp->mgc1[iter]-10][0][q]==1)
        {
        val+=gcfeld[p->mgc1[iter]-10][0][q];
        ++count;
        }
        
        // 4
        for(q=1;q<=3;++q)
        if(p->gcorig1[pp->mgc1[iter]-10][3][q]==1)
        {
        val+=gcfeld[p->mgc1[iter]-10][3][q];
        ++count;
        }
        
        // 3
        for(q=1;q<=3;++q)
        if(p->gcorig1[pp->mgc1[iter]-10][2][q]==1)
        {
        val+=gcfeld[p->mgc1[iter]-10][2][q];
        ++count;
        }
        
        // 2
        for(q=1;q<=3;++q)
        if(p->gcorig1[pp->mgc1[iter]-10][1][q]==1)
        {
        val+=gcfeld[p->mgc1[iter]-10][1][q];
        ++count;
        }
        
        // 5
        for(q=1;q<=3;++q)
        if(p->gcorig1[pp->mgc1[iter]-10][4][q]==1)
        {
        val+=gcfeld[p->mgc1[iter]-10][4][q];
        ++count;
        }
        
        // 6
        for(q=1;q<=3;++q)
        if(p->gcorig1[pp->mgc1[iter]-10][5][q]==1)
        {
        val+=gcfeld[p->mgc1[iter]-10][5][q];
        ++count;
        }

		 
        if(count>0)
        val/=double(count);
		 
		V[iter]=val;	 
	 }
}









