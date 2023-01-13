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

#include"resize.h"
#include"lexer.h"

void resize_class::Dresize(double *&field, int iold, int inew)
{
	double *cache;
	int i;
	
	Darray(cache,inew);
	
	for(i=0;i<MIN(iold,inew);++i)
	cache[i]=field[i];
	
	del_Darray(field,iold);
	
	field=cache;
}

void resize_class::Dresize(double **&field, int iold, int inew, int jold, int jnew)
{
	double **cache;
	int i,j;
	
	Darray(cache,inew,jnew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold,jnew);++j)
	cache[i][j]=field[i][j];
	
	del_Darray(field,iold,jold);
	
	field=cache;
}

void resize_class::Dresize(double ***&field, int iold, int inew, int jold, int jnew, int kold, int knew)
{
	double ***cache;
	int i,j,k;
	
	Darray(cache,inew,jnew,knew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold,jnew);++j)
	for(k=0;k<MIN(kold,knew);++k)
	cache[i][j][k]=field[i][j][k];
	
	del_Darray(field,iold,jold,kold);
	
	field=cache;
}

void resize_class::Dresize(double ****&field, int iold, int inew, int jold, int jnew, int kold, int knew, int lold, int lnew)
{
	double ****cache;
	int i,j,k,l;
	
	Darray(cache,inew,jnew,knew,lnew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold,jnew);++j)
	for(k=0;k<MIN(kold,knew);++k)
	for(l=0;l<MIN(lold,lnew);++l)
	cache[i][j][k][l]=field[i][j][k][l];
	
	del_Darray(field,iold,jold,kold,lold);
	
	field=cache;
}

void resize_class::Iresize(int *&field, int iold, int inew)
{
	int *cache;
	int i;
	
	Iarray(cache,inew);
	
	for(i=0;i<MIN(iold,inew);++i)
	cache[i]=field[i];
	
	del_Iarray(field,iold);
	
	field=cache;
}

void resize_class::Iresize(int **&field, int iold, int inew, int jold, int jnew)
{
	int **cache;
	int i,j;
	
	Iarray(cache,inew,jnew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold,jnew);++j)
	cache[i][j]=field[i][j];
	
	del_Iarray(field,iold,jold);

	field=cache;
}

void resize_class::Iresize(int ***&field, int iold, int inew, int jold, int jnew, int kold, int knew)
{
	int ***cache;
	int i,j,k;
	
	Iarray(cache,inew,jnew,knew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold,jnew);++j)
	for(k=0;k<MIN(kold,knew);++k)
	cache[i][j][k]=field[i][j][k];
	
	del_Iarray(field,iold,jold,kold);
	
	field=cache;
}

void resize_class::Iresize(int ***&field, int iold, int inew, int *jold, int *jnew, int kold, int knew)
{
	int ***cache;
	int i,j,k;
	
	Iarray(cache,inew,jnew,knew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold[i],jnew[i]);++j)
	for(k=0;k<MIN(kold,knew);++k)
	cache[i][j][k]=field[i][j][k];
	
	del_Iarray(field,iold,jold,kold);
	
	field=cache;
}

void resize_class::Iresize(int ****&field, int iold, int inew, int jold, int jnew, int kold, int knew, int lold, int lnew)
{
	int ****cache;
	int i,j,k,l;
	
	Iarray(cache,inew,jnew,knew,lnew);
	
	for(i=0;i<MIN(iold,inew);++i)
	for(j=0;j<MIN(jold,jnew);++j)
	for(k=0;k<MIN(kold,knew);++k)
	for(l=0;l<MIN(lold,lnew);++l)
	cache[i][j][k][l]=field[i][j][k][l];
	
	del_Iarray(field,iold,jold,kold,lold);
	
	field=cache;
}
