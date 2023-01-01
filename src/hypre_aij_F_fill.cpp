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
    pgc->flagx7(p,rownum7);

	n=0;
	FLOOP
	{
	count=0;
	val[count] = M.p[n];	
	col[count] = rownum7[FIJK];
	rownum = rownum7[FIJK];
	++count;

	
    if(p->flag7[FIm1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum7[FIm1JK];
	++count;
	}
    
    if(p->flag7[FIp1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum7[FIp1JK];
	++count;
	}
    
    if(p->flag7[FIJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum7[FIJm1K];
	++count;
	}
    
    if(p->flag7[FIJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum7[FIJp1K];
	++count;
	}
    
    if(p->flag7[FIJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum7[FIJKm1];
	++count;
	}
    
    if(p->flag7[FIJKp1]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum7[FIJKp1];
	++count;
	}
	
	HYPRE_IJMatrixSetValues(A, 1, &count, &rownum, col, val);
	
	++n;
	}
	
	
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	
    // vec
	n=0;
	FLOOP
	{
		xvec[n] = f[FIJK];
		rows[n] = rownum7[FIJK];
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

void hypre_aij::fill_matrix_F_7p_v2(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    fieldint4 rownum4(p);
    
    pgc->rownum4_update(p,rownum4);
    pgc->facenbx(p,rownum4,p->range_row4);

	n=0;
	LOOP
	{
	count=0;
	val[count] = M.p[n];	
	col[count] = rownum4(i,j,k);
	rownum = rownum4(i,j,k);
	++count;

	
    if(p->flag7[FIm1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum4(i-1,j,k);
	++count;
	}
    
    if(p->flag7[FIp1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum4(i+1,j,k);
	++count;
	}
    
    if(p->flag7[FIJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum4(i,j-1,k);
	++count;
	}
    
    if(p->flag7[FIJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum4(i,j+1,k);
	++count;
	}
    
    if(p->flag7[FIJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum4(i,j,k-1);
	++count;
	}
    
    if(p->flag7[FIJKp1]>0 && p->flag7[FIJKp2]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum4(i,j,k+1);
	++count;
	}
	
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
		rows[n] = rownum4(i,j,k);
	++n;
	}

	HYPRE_IJVectorSetValues(b, p->N4_row, rows, rhsvec.V);
	HYPRE_IJVectorSetValues(x, p->N4_row, rows, xvec);
	
	HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
    
}


void hypre_aij::fill_matrix_F_13p(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    int* rownum7;
    p->Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
    
    pgc->rownum7_update(p,rownum7);
    pgc->flagx7(p,rownum7);

    
	n=0;
	FLOOP
	{
	count=0;
	
	val[count] = M.p[n];
	col[count] = rownum7[FIJK];
	rownum = rownum7[FIJK];
	++count;
	
	
    if(p->flag7[FIm1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum7[FIm1JK];
	++count;
	}
    
    if(p->flag7[FIp1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum7[FIp1JK];
	++count;
	}
    
    if(p->flag7[FIJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum7[FIJm1K];
	++count;
	}
    
    if(p->flag7[FIJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum7[FIJp1K];
	++count;
	}
    
    if(p->flag7[FIJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum7[FIJKm1];
	++count;
	}
    
    if(p->flag7[FIJKp1]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum7[FIJKp1];
	++count;
	}
    
    
    // -- 
    if(p->flag7[FIm2JK]>0)
	{
	val[count] = M.ss[n];
	col[count] = rownum7[FIm2JK];
	++count;
	}
    
    if(p->flag7[FIp2JK]>0)
	{
	val[count] = M.nn[n];
	col[count] = rownum7[FIp2JK];
	++count;
	}
    
    if(p->flag7[FIJm2K]>0)
	{
	val[count] = M.ee[n];
	col[count] = rownum7[FIJm2K];
	++count;
	}
    
    if(p->flag7[FIJp2K]>0)
	{
	val[count] = M.ww[n];
	col[count] = rownum7[FIJp2K];
	++count;
	}
    
    if(p->flag7[FIJKm2]>0)
	{
	val[count] = M.bb[n];
	col[count] = rownum7[FIJKm2];
	++count;
	}
    
    if(p->flag7[FIJKp2]>0)
	{
	val[count] = M.tt[n];
	col[count] = rownum7[FIJKp2];
	++count;
	}

	
	HYPRE_IJMatrixSetValues(A, 1, &count, &rownum, col, val);
	
	++n;
	}
	
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	
    // vec
	n=0;
	FLOOP
	{
		xvec[n] = f[FIJK];
		rows[n] = rownum7[FIJK];
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

void hypre_aij::fill_matrix_F_13p_v2(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    fieldint4 rownum4(p);
    
    pgc->rownum4_update(p,rownum4);
    pgc->gcparaxint(p,rownum4,4);

    
	n=0;
	LOOP
	{
	count=0;
	
	val[count] = M.p[n];
	col[count] = rownum4(i,j,k);
	rownum = rownum4(i,j,k);
	++count;
	
	
    if(p->flag7[FIm1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum4(i-1,j,k);
	++count;
	}
    
    if(p->flag7[FIp1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum4(i+1,j,k);
	++count;
	}
    
    if(p->flag7[FIJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum4(i,j-1,k);
	++count;
	}
    
    if(p->flag7[FIJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum4(i,j+1,k);
	++count;
	}
    
    if(p->flag7[FIJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum4(i,j,k-1);
	++count;
	}
    
    if(p->flag7[FIJKp1]>0 && p->flag7[FIJKp2]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum4(i,j,k+1);
	++count;
	}
    
    
    // -- 
    if(p->flag7[FIm2JK]>0)
	{
	val[count] = M.ss[n];
	col[count] = rownum4(i-2,j,k);
	++count;
	}
    
    if(p->flag7[FIp2JK]>0)
	{
	val[count] = M.nn[n];
	col[count] = rownum4(i+2,j,k);
	++count;
	}
    
    if(p->flag7[FIJm2K]>0)
	{
	val[count] = M.ee[n];
	col[count] = rownum4(i,j-2,k);
	++count;
	}
    
    if(p->flag7[FIJp2K]>0)
	{
	val[count] = M.ww[n];
	col[count] = rownum4(i,j+2,k);
	++count;
	}
    
    if(p->flag7[FIJKm2]>0)
	{
	val[count] = M.bb[n];
	col[count] = rownum4(i,j,k-2);
	++count;
	}
    
    if(p->flag7[FIJKp2]>0 && p->flag7[FIJKp3]>0)
	{
	val[count] = M.tt[n];
	col[count] = rownum4(i,j,k+2);
	++count;
	}
    
    //

	
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
		rows[n] = rownum4(i,j,k);
	++n;
	}

	HYPRE_IJVectorSetValues(b, p->N4_row, rows, rhsvec.V);
	HYPRE_IJVectorSetValues(x, p->N4_row, rows, xvec);
	
	HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
}

void hypre_aij::fill_matrix_F_19p(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    int* rownum7;
    p->Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
    
    pgc->rownum7_update(p,rownum7);
    pgc->flagx7(p,rownum7);

    
	n=0;
	FLOOP
	{
	count=0;
	
	val[count] = M.p[n];
	col[count] = rownum7[FIJK];
	rownum = rownum7[FIJK];
	++count;
	
	
    if(p->flag7[FIm1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum7[FIm1JK];
	++count;
	}
    
    if(p->flag7[FIp1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum7[FIp1JK];
	++count;
	}
    
    if(p->flag7[FIJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum7[FIJm1K];
	++count;
	}
    
    if(p->flag7[FIJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum7[FIJp1K];
	++count;
	}
    
    if(p->flag7[FIJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum7[FIJKm1];
	++count;
	}
    
    if(p->flag7[FIJKp1]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum7[FIJKp1];
	++count;
	}
    
    
    // -- 
    if(p->flag7[FIm2JK]>0)
	{
	val[count] = M.ss[n];
	col[count] = rownum7[FIm2JK];
	++count;
	}
    
    if(p->flag7[FIp2JK]>0)
	{
	val[count] = M.nn[n];
	col[count] = rownum7[FIp2JK];
	++count;
	}
    
    if(p->flag7[FIJm2K]>0)
	{
	val[count] = M.ee[n];
	col[count] = rownum7[FIJm2K];
	++count;
	}
    
    if(p->flag7[FIJp2K]>0)
	{
	val[count] = M.ww[n];
	col[count] = rownum7[FIJp2K];
	++count;
	}
    
    if(p->flag7[FIJKm2]>0)
	{
	val[count] = M.bb[n];
	col[count] = rownum7[FIJKm2];
	++count;
	}
    
    if(p->flag7[FIJKp2]>0)
	{
	val[count] = M.tt[n];
	col[count] = rownum7[FIJKp2];
	++count;
	}
    
    
    // -- 
    if(p->flag7[FIm3JK]>0)
	{
	val[count] = M.sss[n];
	col[count] = rownum7[FIm3JK];
	++count;
	}
    
    if(p->flag7[FIp3JK]>0)
	{
	val[count] = M.nnn[n];
	col[count] = rownum7[FIp3JK];
	++count;
	}
    
    if(p->flag7[FIJm3K]>0)
	{
	val[count] = M.eee[n];
	col[count] = rownum7[FIJm3K];
	++count;
	}
    
    if(p->flag7[FIJp3K]>0)
	{
	val[count] = M.www[n];
	col[count] = rownum7[FIJp3K];
	++count;
	}
    
    if(p->flag7[FIJKm3]>0)
	{
	val[count] = M.bbb[n];
	col[count] = rownum7[FIJKm3];
	++count;
	}
    
    if(p->flag7[FIJKp3]>0)
	{
	val[count] = M.ttt[n];
	col[count] = rownum7[FIJKp3];
	++count;
	}
	
	HYPRE_IJMatrixSetValues(A, 1, &count, &rownum, col, val);
	
	++n;
	}
	
	
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	
    // vec
	n=0;
	FLOOP
	{
		xvec[n] = f[FIJK];
		rows[n] = rownum7[FIJK];
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

void hypre_aij::fill_matrix_F_19p_v2(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    fieldint4 rownum4(p);
    
    pgc->rownum4_update(p,rownum4);
    pgc->gcparaxint(p,rownum4,4);

    
	n=0;
	LOOP
	{
	count=0;
	
	val[count] = M.p[n];
	col[count] = rownum4(i,j,k);
	rownum = rownum4(i,j,k);
	++count;
	
	
    if(p->flag7[FIm1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum4(i-1,j,k);
	++count;
	}
    
    if(p->flag7[FIp1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum4(i+1,j,k);
	++count;
	}
    
    if(p->flag7[FIJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum4(i,j-1,k);
	++count;
	}
    
    if(p->flag7[FIJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum4(i,j+1,k);
	++count;
	}
    
    if(p->flag7[FIJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum4(i,j,k-1);
	++count;
	}
    
    if(p->flag7[FIJKp1]>0 && p->flag7[FIJKp2]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum4(i,j,k+1);
	++count;
	}
    
    
    // -- 
    if(p->flag7[FIm2JK]>0)
	{
	val[count] = M.ss[n];
	col[count] = rownum4(i-2,j,k);
	++count;
	}
    
    if(p->flag7[FIp2JK]>0)
	{
	val[count] = M.nn[n];
	col[count] = rownum4(i+2,j,k);
	++count;
	}
    
    if(p->flag7[FIJm2K]>0)
	{
	val[count] = M.ee[n];
	col[count] = rownum4(i,j-2,k);
	++count;
	}
    
    if(p->flag7[FIJp2K]>0)
	{
	val[count] = M.ww[n];
	col[count] = rownum4(i,j+2,k);
	++count;
	}
    
    if(p->flag7[FIJKm2]>0)
	{
	val[count] = M.bb[n];
	col[count] = rownum4(i,j,k-2);
	++count;
	}
    
    if(p->flag7[FIJKp2]>0 && p->flag7[FIJKp3]>0)
	{
	val[count] = M.tt[n];
	col[count] = rownum4(i,j,k+2);
	++count;
	}
    
    // -- 
    if(p->flag7[FIm3JK]>0)
	{
	val[count] = M.sss[n];
	col[count] = rownum4(i-3,j,k);
	++count;
	}
    
    if(p->flag7[FIp3JK]>0)
	{
	val[count] = M.nnn[n];
	col[count] = rownum4(i+3,j,k);
	++count;
	}
    
    if(p->flag7[FIJm3K]>0)
	{
	val[count] = M.eee[n];
	col[count] = rownum4(i,j-3,k);
	++count;
	}
    
    if(p->flag7[FIJp3K]>0)
	{
	val[count] = M.www[n];
	col[count] = rownum4(i,j+3,k);
	++count;
	}
    
    if(p->flag7[FIJKm3]>0)
	{
	val[count] = M.bbb[n];
	col[count] = rownum4(i,j,k-3);
	++count;
	}
    
    if(p->flag7[FIJKp3]>0 && p->flag7[FIJKp4]>0)
	{
	val[count] = M.ttt[n];
	col[count] = rownum4(i,j,k+3);
	++count;
	}

	
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
		rows[n] = rownum4(i,j,k);
	++n;
	}

	HYPRE_IJVectorSetValues(b, p->N4_row, rows, rhsvec.V);
	HYPRE_IJVectorSetValues(x, p->N4_row, rows, xvec);
	
	HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
}

#endif
