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

#include"hypre_struct2D.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sliceint4.h"

void hypre_struct2D::fill_matrix(lexer* p, ghostcell* pgc, matrix2D &M, slice &f, vec2D &rhsvec)
{
    sliceint4 cval4(p);
    
    count=0;
    SLICELOOP4
    {
    cval4(i,j)=count;
    ++count;
    }
    
    nentries=5;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

    count=0;
    JILOOP
    {
		PSLICECHECK4
		{
		n=cval4(i,j);
        
		values[count]=M.p[n];
		++count;
		
		values[count]=M.s[n];
		++count;
		
		values[count]=M.n[n];
		++count;
		
		values[count]=M.e[n];
		++count;
		
		values[count]=M.w[n];
		++count;
		}     
		
        SSLICECHECK4
		{
		values[count]=1.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		}    
    }
	
    HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, stencil_indices, values);
    HYPRE_StructMatrixAssemble(A);
    
    
    // vec
    count=0;
	JILOOP
	{
		PSLICECHECK4
		values[count] = f(i,j);
		
		SSLICECHECK4
		values[count] = 0.0;
	
    ++count;
    }

    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
    HYPRE_StructVectorAssemble(x);
    
    
    count=0; 
	JILOOP
	{
		PSLICECHECK4
		{
		n=cval4(i,j);
		values[count] = rhsvec.V[n];
		}
		
		SSLICECHECK4
		values[count] = 0.0;

    ++count;
    }
    
    HYPRE_StructVectorSetBoxValues(rhs, ilower, iupper, values);
    HYPRE_StructVectorAssemble(rhs);
}

void hypre_struct2D::fillbackvec(lexer *p, slice &f, vec2D &xvec, int var)
{
	HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);
	
        count=0;
        JILOOP
        {
		PSLICECHECK4
        f(i,j)=values[count];
		
        ++count;
        }
}

#endif
