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

#include"potentialfile_in.h"
#include"lexer.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>


void potentialfile_in::header_read(lexer *p, ghostcell *pgc)
{
    // Open File
    sprintf(name,"./REEF3D_PotentialFile/REEF3D-potentialheader-%i.r3d",p->I240);
	
    
    // count entries
	headerfile.open(name, ios::binary);
    
    entrycount=0;
    
    while(!headerfile.eof())
	{
    headerfile.read((char*)&iin, sizeof (int));
    headerfile.read((char*)&ffn, sizeof (float));
    
    ++entrycount;
    }
    headerfile.close();

    //alloocate arrays
    p->Iarray(iter,entrycount);
    p->Darray(xloc,entrycount);
    
    // read entries
	headerfile.open(name, ios::binary);
    
    n=0;
    while(!headerfile.eof())
	{
    headerfile.read((char*)&iin, sizeof (int));
    iter[n]=iin;
    headerfile.read((char*)&ffn, sizeof (float));
    xloc[n]=ffn;
    
    ++n;
    }
    headerfile.close();



    

}
