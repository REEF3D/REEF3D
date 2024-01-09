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

#include"hypre_struct_fnpf.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fieldint4.h"
#include"matrix_diag.h"

void hypre_struct_fnpf::fill_matrix8_2Dvert(lexer* p, ghostcell* pgc, double *f, double *rhs, double *M)
{
    nentries=9;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

    // M
    HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, stencil_indices, M);
    HYPRE_StructMatrixAssemble(A);
    
    // x
    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, f);
    HYPRE_StructVectorAssemble(x);

    // rhs
    HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs);
    HYPRE_StructVectorAssemble(b);
    
}

#endif
