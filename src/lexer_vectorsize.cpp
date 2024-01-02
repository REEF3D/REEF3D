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
#include"ghostcell.h"

void lexer::vecsize(ghostcell *pgc)
{
    int n;
    gcbextra=10;
	int gcbnum=0;
    int safetymargin=0;
    int solid_gcb_est_max, topo_gcb_est_max, gcextra_max;
    
    int gcbextra0;

    solid_gcb_est_max = pgc->globalimax(solid_gcb_est);
    topo_gcb_est_max = pgc->globalimax(topo_gcb_est);
    
    gcextra_max = pgc->globalimax(gcextra4);

	gcb_sediment_est = gcb4_count*margin;	
    
    gcb_floating_est=0;
    
    
// gcbextra
    gcbextra=gcextra_max*margin;
    
    gcbextra0=gcbextra;
    
    safetymargin = 0.2*double(solid_gcbextra_est+topo_gcbextra_est+tot_gcbextra_est) + 100.0;
    
    

    // solid and topo
	if(S10>0 || G1>0)
    gcbextra+=MAX(MAX(solid_gcbextra_est,topo_gcbextra_est),tot_gcbextra_est) + int(safetymargin);
    
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
    

    veclength = cellnum + gcbnum*margin + gcpara_sum*4  + gcbextra;
    
    //cout<<mpirank<<" CELLNUM: "<<cellnum<<" veclength: "<<veclength<<" gcextra0: "<<gcbextra0<<" gcextra: "<<gcbextra<<" tot_gcbextra_est: "<<tot_gcbextra_est<<endl;
    
    //gcbextra=gcbextra0;
    
    C4_size=C4a_size=C6_size=M_size=veclength;
    
    
    
    
    // 2D
     slicenum= 0;
    for(i=0; i<knox; ++i)
    for(j=0; j<knoy; ++j)
	if(flagslice4[(i-imin)*jmax + (j-jmin)]>0)
    ++slicenum;
    
    gcpara_sum = gcslpara1_count + gcslpara2_count + gcslpara3_count + gcslpara4_count
              + gcslparaco1_count + gcslparaco2_count + gcslparaco3_count + gcslparaco4_count;
    
    vec2Dlength = slicenum + gcbnum*3  + gcpara_sum*4;    
    
    C1_2D_size=C2_2D_size=C4_2D_size=M_2D_size=vec2Dlength;
}

void lexer::gcbextra_est(ghostcell *pgc)
{
    
    
    
    
    
}
