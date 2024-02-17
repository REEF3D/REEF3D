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

#include"potentialfile_out.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void potentialfile_out::header_file_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->mpirank==0)
    for(n=0;n<p->P240;++n)
    {
    // open file
    sprintf(headername,"./REEF3D_FlowFile/REEF3D-potentialheader.r3d");
		
    // openfile
    headerout.open(headername, ios::binary);
    }
    
     // header
    if(p->mpirank==0)
    {
        iin=p->P240;
        headerout.write((char*)&iin, sizeof (int));
        
        for(n=0;n<p->P240;++n)
        {  
            iin=n;
            headerout.write((char*)&iin, sizeof (int));
            ffn=float(p->P240_x[n]);
            headerout.write((char*)&ffn, sizeof (float));

        headerout.close();
        }
    }
}

void potentialfile_out::header_file(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{   
}

