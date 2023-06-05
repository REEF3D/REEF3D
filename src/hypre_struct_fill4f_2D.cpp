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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_struct::fill_matrix4f_2Dvert(lexer* p, ghostcell* pgc, field &f, vec &rhs, matrix_diag &M)
{
    count=0;
    FLUIDLOOP
    {
    CVAL4[IJK]=count;
    ++count;
    }
    
    nentries=5;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

    count=0;
    KJILOOP
    {
		PFLUIDCHECK
		{
		n=CVAL4[IJK];
        
		values[count]=M.p[n];
		++count;
		
		values[count]=M.s[n];
		++count;
		
		values[count]=M.n[n];
		++count;

		values[count]=M.b[n];
		++count;
		
		values[count]=M.t[n];
		++count; 
		}     
		
		SFLUIDCHECK
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
	KJILOOP
	{
		PFLUIDCHECK
		values[count] = f(i,j,k);
		
		SFLUIDCHECK
		values[count] = 0.0;
	
    ++count;
    }

    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
    HYPRE_StructVectorAssemble(x);
    
    
    count=0; 
	KJILOOP
	{
		PFLUIDCHECK
		{
		n=CVAL4[IJK];
		values[count] = rhs.V[n];
		}
		
		SFLUIDCHECK
		values[count] = 0.0;

    ++count;
    }
    
    HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);
    HYPRE_StructVectorAssemble(b);
}


#endif
