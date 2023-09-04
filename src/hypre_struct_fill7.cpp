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
#include"ghostcell.h"
#include"matrix_diag.h"

void hypre_struct::fill_matrix7(lexer* p, ghostcell* pgc, double *f, vec &rhs, matrix_diag &M)
{

    count=0;
    LOOP
    {
    CVAL4[IJK]=count;
    ++count;
    }
    
    nentries=7;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

    count=0;
    KJILOOP
    {
		PCHECK
		{
		n=CVAL4[IJK];
        
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
		
		values[count]=M.b[n];
		++count;
		
		values[count]=M.t[n];
		++count; 
		}     
		
		SCHECK
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
		PCHECK
		values[count] = f[IJK];
		
		SCHECK
		values[count] = 0.0;
	
    ++count;
    }

    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
    HYPRE_StructVectorAssemble(x);
    
    
    count=0; 
	KJILOOP
	{
		FPCHECK
		{
		n=CVAL4[IJK];
		values[count] = rhs.V[n];
		}
		
		FSCHECK
		values[count] = 0.0;

    ++count;
    }
    
    HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);
    HYPRE_StructVectorAssemble(b);

}

void hypre_struct::fillbackvec7(lexer *p, double *f, int var)
{
	HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);
	
        count=0;
        KJILOOP
        {
        PCHECK
        f[IJK]=values[count];
		
        ++count;
        }
}

#endif
