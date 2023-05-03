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

#include"hypre_sstruct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_sstruct::fill_matrixM(lexer* p, ghostcell* pgc, double *f, double *rhs, double *M)
{  
    nentries=15;
    
    for (j=0; j< nentries; j++)
    stencil_indices[j] = j;

	// M
    HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, variable, nentries, stencil_indices, M);
    HYPRE_SStructMatrixAssemble(A);
    
    // x
    HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, variable, f);
    HYPRE_SStructVectorAssemble(x);
    
    // b
    HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, variable, rhs);
    HYPRE_SStructVectorAssemble(b);
}

void hypre_sstruct::fill_matrixM_2Dvert(lexer* p, ghostcell* pgc, double *f, double *rhs, double *M)
{
    nentries=9;
    
    for (j=0; j< nentries; j++)
    stencil_indices[j] = j;

    // M
    HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, variable, nentries, stencil_indices, M);
    HYPRE_SStructMatrixAssemble(A);
    
    // x
    HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, variable, f);
    HYPRE_SStructVectorAssemble(x);

    // rhs
    HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, variable, rhs);
    HYPRE_SStructVectorAssemble(b);
}

void hypre_sstruct::fillbackvecM(lexer *p, double *f)
{
	HYPRE_SStructVectorGetBoxValues(x, part, ilower, iupper, variable, f);
}

#endif
