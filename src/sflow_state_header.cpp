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

#include"sflow_state.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void sflow_state::write_header(lexer *p, fdm2D *b, ghostcell *pgc)
{
    ofstream headout;
    
    
    // file name
    filename_header(p,b,pgc);
    
    // open file
	headout.open(name, ios::binary);
    
    // ini write    
    iin=p->origin_i-is_global;
    
    if(is_flag==1)
    iin=0;
    
    headout.write((char*)&iin, sizeof (int));
    
    // ---
    iin=p->origin_j-js_global;
    
    if(js_flag==1)
    iin=0;
    
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->origin_k;
    headout.write((char*)&iin, sizeof (int));


    ddn=p->originx;
    headout.write((char*)&ddn, sizeof (double));
    
    ddn=p->originy;
    headout.write((char*)&ddn, sizeof (double));
    
    ddn=p->originz;
    headout.write((char*)&ddn, sizeof (double));
  
    
    iin=ie-is;
    headout.write((char*)&iin, sizeof (int));
    
    iin=je-js;
    headout.write((char*)&iin, sizeof (int));
    
    iin=1; //knoz
    headout.write((char*)&iin, sizeof (int));
    
    
    iin=p->nb1;
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->nb2;
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->nb3;
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->nb4;
    headout.write((char*)&iin, sizeof (int));
    
    
    //
    for(i=is;i<ie;++i)
    {
    ddn=p->XP[IP];
    headout.write((char*)&ddn, sizeof (double));
    } 
    
    for(j=js;j<je;++j)
    {
    ddn=p->YP[JP];
    headout.write((char*)&ddn, sizeof (double));
    } 
    
    // ZN
    ddn=0.0;
    headout.write((char*)&ddn, sizeof (double));

    
    for(i=is;i<ie;++i)
    for(j=js;j<je;++j)
    PSLICECHECK4
    {
    ffn=float(b->bed(i,j));
    headout.write((char*)&ffn, sizeof (float));
    } 
    
    headout.close();
}


