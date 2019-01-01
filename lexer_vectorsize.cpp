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

void lexer::vecsize()
{
    /*
    int *flag;
    
    
    Iarray(flag,imax*jmax*kmax);
    
    
    for(i=0;i<imax*jmax*kmax;++i)
    flag[i] = flag4[i];
    
    
    LOOP
    if(flag_solid[IJK] < 0 && G39==1 )
    flag[IJK] = -1;
    
    LOOP
    if(flag_solid[IJK] < 0 && G39==1 )
    flag[IJK] = -1;
    
    */
    int n;
    gcbextra=0;
	int gcbnum=0;


	gcb_sediment_est = gcb4_count*3;	
	gcb_floating_est = gcb4_count;
    
    gcbextra=gcextra4*3;
    
    // solid and topo
	if(S10>0 || G1>0)
	gcbextra+=(solid_gcb_est*4+geotopo_gcb_est);
    
    // floating 
	if(X10>0)
	gcbextra+=(gcb_floating_est);
    
    // topo for sediment
    if(S10>0 )
	gcbextra+=gcb_sediment_est + knoy*knox;
    
    // extra allocate
    gcbextra += int(cellnum/75);
	

    stencil=7;
	
	gcbnum = MAX(gcbnum,gcb1_count);
	gcbnum = MAX(gcbnum,gcb2_count);
	gcbnum = MAX(gcbnum,gcb3_count);
	gcbnum = MAX(gcbnum,gcb4_count);
	
	gcbnum+=500;
    

    //cout<<mpirank<<" CELLNUM: "<<cellnum<<endl;
    
    veclength = cellnum + gcbnum*3 + gcpara_sum*4  + gcbextra;
    
    C1_size=C2_size=C3_size=C4_size=C4a_size=C6_size=M_size=veclength;
}

