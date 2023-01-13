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

#include"hypre_struct_fnpf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"ghostcell.h"
#include"fieldint4.h"
#include"matrix_diag.h"

void hypre_struct_fnpf::fill_matrix8(lexer* p, ghostcell* pgc, double *f, double *rhs, double *M)
{    
    nentries=15;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

   
	// M
    HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, stencil_indices, M);
    HYPRE_StructMatrixAssemble(A);
    
    
    // x
    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, f);
    HYPRE_StructVectorAssemble(x);
    
    // b
    HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs);
    HYPRE_StructVectorAssemble(b);
    
}

void hypre_struct_fnpf::fillbackvec8(lexer *p, double *f, double *rhs, double *M)
{
	HYPRE_StructVectorGetBoxValues(x, ilower, iupper, f);
}

#endif
