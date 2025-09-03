/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"hypre_aij.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void hypre_aij::fill_matrix_F_7p(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    int* rownum7;
    p->Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
    
    pgc->rownum7_update(p,rownum7);
    pgc->flagx(p,rownum7);

	n=0;
	LOOP
	{
	count=0;
	val[count] = M.p[n];	
	col[count] = rownum7[IJK];
	rownum = rownum7[IJK];
	++count;

	
	val[count] = M.s[n];
	col[count] = rownum7[Im1JK];
	++count;
    
	val[count] = M.n[n];
	col[count] = rownum7[Ip1JK];
	++count;

	val[count] = M.e[n];
	col[count] = rownum7[IJm1K];
	++count;
 
	val[count] = M.w[n];
	col[count] = rownum7[IJp1K];
	++count;

	val[count] = M.b[n];
	col[count] = rownum7[IJKm1];
	++count;

	val[count] = M.t[n];
	col[count] = rownum7[IJKp1];
	++count;

	
	HYPRE_IJMatrixSetValues(A, 1, &count, &rownum, col, val);
	
	++n;
	}
	
	
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	
    // vec
	n=0;
	LOOP
	{
		xvec[n] = f[FIJK];
		rows[n] = rownum7[IJK];
	++n;
	}

	HYPRE_IJVectorSetValues(b, p->N7_row, rows, rhsvec.V);
	HYPRE_IJVectorSetValues(x, p->N7_row, rows, xvec);
	
	HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
    
    p->del_Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
}

#endif
