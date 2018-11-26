/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"hypre_aij.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void hypre_aij::fill_matrix_F(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
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
	
	
    if(p->flag7[Im1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum7[FIm1JK];
	++count;
	}
    
    if(p->flag7[Ip1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum7[FIp1JK];
	++count;
	}
    
    if(p->flag7[IJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum7[FIJm1K];
	++count;
	}
    
    if(p->flag7[IJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum7[FIJp1K];
	++count;
	}
    
    if(p->flag7[IJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum7[FIJKm1];
	++count;
	}
    
    if(p->flag7[IJKp1]>0)
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


void hypre_aij::fill_matrix_F_13p(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    int* rownum7;
    p->Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
    
    pgc->rownum7_update(p,rownum7);
    pgc->flagx7(p,rownum7);
    
    HYPRE_IJMatrixInitialize(A);
    HYPRE_IJVectorInitialize(x);
    HYPRE_IJVectorInitialize(b);
    
	n=0;
	FLOOP
	{
	count=0;
	
	val[count] = M.p[n];
	col[count] = rownum7[FIJK];
	rownum = rownum7[FIJK];
	++count;
	
	
    if(p->flag7[Im1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum7[FIm1JK];
	++count;
	}
    
    if(p->flag7[Ip1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum7[FIp1JK];
	++count;
	}
    
    if(p->flag7[IJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum7[FIJm1K];
	++count;
	}
    
    if(p->flag7[IJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum7[FIJp1K];
	++count;
	}
    
    if(p->flag7[IJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum7[FIJKm1];
	++count;
	}
    
    if(p->flag7[IJKp1]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum7[FIJKp1];
	++count;
	}
    
    
    // -- 
    if(p->flag7[Im2JK]>0)
	{
	val[count] = M.ss[n];
	col[count] = rownum7[FIm2JK];
	++count;
	}
    
    if(p->flag7[Ip2JK]>0)
	{
	val[count] = M.nn[n];
	col[count] = rownum7[FIp2JK];
	++count;
	}
    
    if(p->flag7[IJm2K]>0)
	{
	val[count] = M.ee[n];
	col[count] = rownum7[FIJm2K];
	++count;
	}
    
    if(p->flag7[IJp2K]>0)
	{
	val[count] = M.ww[n];
	col[count] = rownum7[FIJp2K];
	++count;
	}
    
    if(p->flag7[IJKm2]>0)
	{
	val[count] = M.bb[n];
	col[count] = rownum7[FIJKm2];
	++count;
	}
    
    if(p->flag7[IJKp2]>0)
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

void hypre_aij::fill_matrix_F_19p(lexer* p, ghostcell* pgc, matrix_diag &M, double *f, double *xvec, vec &rhsvec)
{
    int* rownum7;
    p->Iarray(rownum7,p->imax*p->jmax*(p->kmax+2));
    
    pgc->rownum7_update(p,rownum7);
    pgc->flagx7(p,rownum7);
    
    HYPRE_IJMatrixInitialize(A);
    HYPRE_IJVectorInitialize(x);
    HYPRE_IJVectorInitialize(b);
    
	n=0;
	FLOOP
	{
	count=0;
	
	val[count] = M.p[n];
	col[count] = rownum7[FIJK];
	rownum = rownum7[FIJK];
	++count;
	
	
    if(p->flag7[Im1JK]>0)
	{
	val[count] = M.s[n];
	col[count] = rownum7[FIm1JK];
	++count;
	}
    
    if(p->flag7[Ip1JK]>0)
	{
	val[count] = M.n[n];
	col[count] = rownum7[FIp1JK];
	++count;
	}
    
    if(p->flag7[IJm1K]>0)
	{
	val[count] = M.e[n];
	col[count] = rownum7[FIJm1K];
	++count;
	}
    
    if(p->flag7[IJp1K]>0)
	{
	val[count] = M.w[n];
	col[count] = rownum7[FIJp1K];
	++count;
	}
    
    if(p->flag7[IJKm1]>0)
	{
	val[count] = M.b[n];
	col[count] = rownum7[FIJKm1];
	++count;
	}
    
    if(p->flag7[IJKp1]>0)
	{
	val[count] = M.t[n];
	col[count] = rownum7[FIJKp1];
	++count;
	}
    
    
    // -- 
    if(p->flag7[Im2JK]>0)
	{
	val[count] = M.ss[n];
	col[count] = rownum7[FIm2JK];
	++count;
	}
    
    if(p->flag7[Ip2JK]>0)
	{
	val[count] = M.nn[n];
	col[count] = rownum7[FIp2JK];
	++count;
	}
    
    if(p->flag7[IJm2K]>0)
	{
	val[count] = M.ee[n];
	col[count] = rownum7[FIJm2K];
	++count;
	}
    
    if(p->flag7[IJp2K]>0)
	{
	val[count] = M.ww[n];
	col[count] = rownum7[FIJp2K];
	++count;
	}
    
    if(p->flag7[IJKm2]>0)
	{
	val[count] = M.bb[n];
	col[count] = rownum7[FIJKm2];
	++count;
	}
    
    if(p->flag7[IJKp2]>0)
	{
	val[count] = M.tt[n];
	col[count] = rownum7[FIJKp2];
	++count;
	}
    
    
    // -- 
    if(p->flag7[Im3JK]>0)
	{
	val[count] = M.sss[n];
	col[count] = rownum7[FIm3JK];
	++count;
	}
    
    if(p->flag7[Ip3JK]>0)
	{
	val[count] = M.nnn[n];
	col[count] = rownum7[FIp3JK];
	++count;
	}
    
    if(p->flag7[IJm3K]>0)
	{
	val[count] = M.eee[n];
	col[count] = rownum7[FIJm3K];
	++count;
	}
    
    if(p->flag7[IJp3K]>0)
	{
	val[count] = M.www[n];
	col[count] = rownum7[FIJp3K];
	++count;
	}
    
    if(p->flag7[IJKm3]>0)
	{
	val[count] = M.bbb[n];
	col[count] = rownum7[FIJKm3];
	++count;
	}
    
    if(p->flag7[IJKp3]>0)
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

#endif