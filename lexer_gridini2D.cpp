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

#include"lexer.h"

void lexer::flagini2D()
{
    control_calc();

	grid2Dsize();
	vellast();
	
	if(N5==0)
	i_dir=j_dir=k_dir=1;
	
	x_dir=y_dir=z_dir=1.0;
	
	if(i_dir==0)
	x_dir=0.0;
	
	if(j_dir==0)
	y_dir=0.0;
	
	if(k_dir==0)
	z_dir=0.0;
	
	
}

void lexer::gridini2D()
{
    Iarray(sizeS1, 5);
    Iarray(sizeS2, 5);
    Iarray(sizeS4, 5);

    for(int n=0;n<5;++n)
    {
    sizeS1[n]=0;
    sizeS2[n]=0;
    sizeS4[n]=0;
    }	
	
// ---	
	int istart,iend,jstart,jend,kstart,kend,qn;
	int count=0;
	
	for(qn=0;qn<G95;++qn)
    {
        istart = conv((G95_xs[qn]-originx)/dx);
        iend = conv((G95_xe[qn]-originx)/dx);

        jstart = conv((G95_ys[qn]-originy)/dx);
        jend = conv((G95_ye[qn]-originy)/dx);

        kstart = conv((G95_zs[qn]-originz)/dx);
        kend = conv((G95_ze[qn]-originz)/dx);

        for(n=0;n<gcb4_count;++n)
		{
		i=gcb4[n][0];
		j=gcb4[n][1];
		k=gcb4[n][2];
		
			if(i>=istart && i<iend && j>=jstart && j<jend && k>=kstart && k<kend && gcb4[n][3]==5 && (gcb4[n][4]==21||gcb4[n][4]==22))
			{
			++count;
			gcb4[n][4]=2;
			}
		}
    }
	count+=gcout_count;
	
	Iresize(gcout, gcout_count,count,6,6);
}


