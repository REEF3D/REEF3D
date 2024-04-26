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

#include"hypre_sstruct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"ghostcell.h"
#include"fieldint4.h"
#include"matrix_diag.h"

void hypre_sstruct::fill_matrix8(lexer* p, ghostcell* pgc, double *f, vec &rhs, matrix_diag &M)
{
    fieldint4 cval4(p);
    
    count=0;
    FLUIDLOOP
    {
    cval4(i,j,k)=count;
    ++count;
    }
    
    nentries=7;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

    count=0;
    KJILOOP
    {
		PFLUIDCHECK
		{
		n=cval4(i,j,k);
        
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
		
		FSCHECK
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
	
    HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, variable, nentries, stencil_indices, values);
    HYPRE_SStructMatrixAssemble(A);
    
    
    // vec
    count=0;
	KJILOOP
	{
		FPCHECK
		values[count] = f[FIJK];
		
		FSCHECK
		values[count] = 0.0;
	
    ++count;
    }

    HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, variable, values);
    HYPRE_SStructVectorAssemble(x);
    
    
    count=0; 
	KJILOOP
	{
		FPCHECK
		{
		n=cval4(i,j,k);
		values[count] = rhs.V[n];
		}
		
		FSCHECK
		values[count] = 0.0;

    ++count;
    }
    
    HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, variable, values);
    HYPRE_SStructVectorAssemble(b);
    
}

void hypre_sstruct::fillbackvec8(lexer *p, double *f, int var)
{
	HYPRE_SStructVectorGetBoxValues(x, part, ilower, iupper, variable, values);
	
        count=0;
        KJILOOP
        {
		 FPCHECK
        f[FIJK]=values[count];
		
        ++count;
        }
}

#endif
