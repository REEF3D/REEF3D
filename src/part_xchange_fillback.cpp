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

#include"part.h"
#include"lexer.h"
#include"ghostcell.h"

void part::xchange_fillback(lexer *p, ghostcell *pgc, double *F)
{
    index_empty = index_empty0;
    
    // fill recv into F
    for(q=0;q<recvnum[0];++q)
    {
    F[Empty[index_empty]] = recv[0][q]; 
    --index_empty;
    }
    
    for(q=0;q<recvnum[1];++q)
    {
    F[Empty[index_empty]] = recv[1][q]; 
    --index_empty;
    }
    
    for(q=0;q<recvnum[2];++q)
    {
    F[Empty[index_empty]] = recv[2][q]; 
    --index_empty;
    }
    
    for(q=0;q<recvnum[3];++q)
    {
    F[Empty[index_empty]] = recv[3][q]; 
    --index_empty;
    }
    
    for(q=0;q<recvnum[4];++q)
    {
    F[Empty[index_empty]] = recv[4][q]; 
    --index_empty;
    }
    
    for(q=0;q<recvnum[5];++q)
    {
    F[Empty[index_empty]] = recv[5][q]; 
    --index_empty;
    }

}

