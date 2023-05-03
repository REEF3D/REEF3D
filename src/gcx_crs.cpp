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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::rangex(lexer* p, int* range, int cellcount)
{

    int *colnum;
    
    p->Iarray(colnum, p->M10+1);
    offset = cellcount;
    
    for(n=0;n<=p->M10;++n)
    p->colnum[n]=0;
	
    MPI_Allgather(&offset,1,MPI_INT,colnum,1,MPI_INT,mpi_comm);

    range[0]=0;

    for(n=0;n<p->M10;++n)
    range[n+1]=range[n]+colnum[n];
    
    p->del_Iarray(colnum, p->M10+1);

}



