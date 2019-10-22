/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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


void lexer::grid2Dsize()
{
	int gcbnum=0;
    int gcparanum=0;

    
	gcbnum = MAX(gcbnum,gcbsl1_count);
	gcbnum = MAX(gcbnum,gcbsl2_count);
	gcbnum = MAX(gcbnum,gcbsl4_count);
	gcbnum = MAX(gcbnum,gcbsl4a_count);
	
	gcbnum+=100;
    
    
    slicenum= 0;
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
	if(flagslice4[(i-imin)*jmax + (j-jmin)]>0)
    ++slicenum;
    
    gcparanum = gcslpara1_count + gcslpara2_count + gcslpara3_count + gcslpara4_count
              + gcslparaco1_count + gcslparaco2_count + gcslparaco3_count + gcslparaco4_count;
    
    vec2Dlength = 2*slicenum + gcbnum*4  + gcparanum*4;    
    
    C1_2D_size=C2_2D_size=C4_2D_size=M_2D_size=vec2Dlength;
    
    cout<<mpirank<<" SLICENUM: "<<slicenum<<" gcbnum: "<<gcbnum<<" gcparanum: "<<gcparanum<<endl;
    
}