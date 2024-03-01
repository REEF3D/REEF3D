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

#include"cfd_state.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void cfd_state::write_header(lexer *p, fdm *a, ghostcell *pgc)
{
    ofstream headout;
    
    
    // file name
    filename_header(p,a,pgc);
    
    // open file
	headout.open(name, ios::binary);
    
    // ijk origin    
    iin=p->origin_i-is_global;
    
    if(is_flag==1)
    iin=0;
    
    headout.write((char*)&iin, sizeof (int));
    
    
    iin=p->origin_j-js_global;
    
    if(js_flag==1)
    iin=0;
    
    headout.write((char*)&iin, sizeof (int));
    
    
    iin=p->origin_k;
    headout.write((char*)&iin, sizeof (int));

    // xyz origin
    ddn=p->originx;
    headout.write((char*)&ddn, sizeof (double));
    
    ddn=p->originy;
    headout.write((char*)&ddn, sizeof (double));
    
    ddn=p->originz;
    headout.write((char*)&ddn, sizeof (double));
  
    // ijk length
    iin=ie-is;
    headout.write((char*)&iin, sizeof (int));
    
    iin=je-js;
    headout.write((char*)&iin, sizeof (int));
    
    
    iin=p->knoz;
    headout.write((char*)&iin, sizeof (int));
    
    // parallel neihbors
    iin=p->nb1;
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->nb2;
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->nb3;
    headout.write((char*)&iin, sizeof (int));
    
    iin=p->nb4;
    headout.write((char*)&iin, sizeof (int));
    
    
    // grid coordinates
    for(i=is;i<=ie;++i)
    {
    ddn=p->XN[IP];
    headout.write((char*)&ddn, sizeof (double));
    } 
    
    for(j=js;j<=je;++j)
    {
    ddn=p->YN[JP];
    headout.write((char*)&ddn, sizeof (double));
    } 
    
    FKLOOP
    {
    ddn=p->ZN[KP];
    headout.write((char*)&ddn, sizeof (double));
    } 

    
    headout.close();
}


