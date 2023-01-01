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

#include"flowfile_out.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void flowfile_out::header_file_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    for(n=0;n<p->P230;++n)
    {
    // open file
    if(p->P14==0)
    sprintf(headername,"REEF3D-flowheader-%i.r3d",n+1);
			
    if(p->P14==1)
    sprintf(headername,"./REEF3D_FlowFile/REEF3D-flowheader-%i.r3d",n+1);
		
    // openfile
    headerout[n].open(headername, ios::binary);
    }
    
     // header
    if(p->mpirank==0)
    for(n=0;n<p->P230;++n)
    {  
        // xs,xe,ys,ye,zs,ze
        ddn=p->global_xmin;
        headerout[n].write((char*)&ddn, sizeof (double));
        ddn=p->global_xmax;
        headerout[n].write((char*)&ddn, sizeof (double));
        ddn=p->global_ymin;
        headerout[n].write((char*)&ddn, sizeof (double));
        ddn=p->global_ymax;
        headerout[n].write((char*)&ddn, sizeof (double));
        ddn=p->global_zmin;
        headerout[n].write((char*)&ddn, sizeof (double));
        ddn=p->global_zmax;
        headerout[n].write((char*)&ddn, sizeof (double));
        
        // dx
        ddn=p->DXM;
        headerout[n].write((char*)&ddn, sizeof (double));
   
        // Ni,Nj,Nk | Ni*Nj*Nk = elnum_all
        iin=Ni;
        headerout[n].write((char*)&iin, sizeof (int));
        iin=Nj;
        headerout[n].write((char*)&iin, sizeof (int));
        iin=Nk;
        headerout[n].write((char*)&iin, sizeof (int));
        
    headerout[n].close();
    }
}

void flowfile_out::header_file(lexer *p, fdm *a, ghostcell *pgc)
{   
    headerout[n].open(headername, ios::binary | ios::app);
    
    iin=p->count;
    headerout[n].write((char*)&iin, sizeof (int));
        
    ddn = p->simtime;
    headerout[n].write((char*)&ddn, sizeof (double));
    
    headerout[n].close();
}

