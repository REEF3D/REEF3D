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

resize_class::resize_class()
{
}

resize_class::~resize_class()
{
}


void resize_class::Darray(double *& field, int numi)
{
    int n;
	
	if(numi>0)
	field = new double[numi];

	for(n=0; n<numi;++n)
	field[n]=0.0;
}

void resize_class::Darray(double **& field, int numi, int numj)
{
	int n,m;
	
	if(numi>0)
	field=new double*[numi];

	for(n=0;n<numi;++n)
	field[n]=new double[numj];

	for(n=0;n<numi;++n)
	for(m=0;m<numj;++m)
	field[n][m]=0.0;
}

void resize_class::Darray(double ***& field, int numi, int numj, int numk)
{
	int n,m,q;
	
	if(numi>0)
	field = new double**[numi];

	for(n=0; n<numi; ++n)
	field[n] = new double*[numj];

	for(n=0; n<numi; ++n)
	for(m=0; m<numj; ++m)
	field[n][m] = new double[numk];

	for(n=0; n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	field[n][m][q]=0.0;
}

void resize_class::Darray(double ****& field, int numi, int numj, int numk, int numl)
{
	int n,m,l,q;
	
	if(numi>0)
	field = new double***[numi];

	for(n=0; n<numi; ++n)
	field[n] = new double**[numj];

	for(n=0; n<numi; ++n)
	for(m=0; m<numj; ++m)
	field[n][m] = new double*[numk];
	
	for(n=0; n<numi; ++n)
	for(m=0; m<numj; ++m)
	for(q=0; q<numk; ++q)
	field[n][m][q] = new double[numl];

	for(n=0;n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	for(l=0;l<numl;++l)
	field[n][m][q][l]=0.0;
}

void resize_class::Darray(double **& field, int numi, int *numj)
{
	int n,m;
	
	if(numi>0)
	field=new double*[numi];

	for(n=0;n<numi;++n)
	field[n]=new double[numj[n]];

	for(n=0; n<numi;++n)
	for(m=0; m<numj[n];++m)
	field[n][m]=0.0;
}

void resize_class::Iarray(int *& field, int numi)
{
	int n;
	
	if(numi>0)
	field=new int[numi];

	for(n=0; n<numi;++n)
	field[n]=0;
}

void resize_class::Iarray(int **& field, int numi, int numj)
{
	int n,m;
	
	if(numi>0)
	field=new int*[numi];

	for(n=0;n<numi;++n)
	field[n]=new int[numj];

	for(n=0; n<numi;++n)
	for(m=0; m<numj;++m)
	field[n][m]=0;
}

void resize_class::Iarray(int ***& field, int numi, int numj, int numk)
{
	int n,m,q;
	
	if(numi>0)
	field=new int**[numi];
	
	for(n=0; n<numi; ++n)
	field[n] = new int*[numj];

	for(n=0; n<numi; ++n)
	for(m=0; m<numj; ++m)
	field[n][m] = new int[numk];

	for(n=0;n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	field[n][m][q]=0;
}

void resize_class::Iarray(int ****& field, int numi, int numj, int numk, int numl)
{
	int n,m,q,r;
	
	if(numi>0)
	field = new int***[numi];

	for(n=0; n<numi; ++n)
	field[n] = new int**[numj];

	for(n=0; n<numi; ++n)
	for(m=0; m<numj; ++m)
	field[n][m] = new int*[numk];

	for(n=0;n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	field[n][m][q]=new int[numl];
	
	for(n=0; n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	for(r=0;r<numl;++r)
	field[n][m][q][r]=0;
}

void resize_class::Iarray(int *****& field, int numi, int numj, int numk, int numl, int numh)
{
	int n,m,q,r,s;
	
	if(numi>0)
	field = new int****[numi];

	for(n=0; n<numi; ++n)
	field[n] = new int***[numj];

	for(n=0; n<numi; ++n)
	for(m=0; m<numj; ++m)
	field[n][m] = new int**[numk];

	for(n=0; n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	field[n][m][q]=new int*[numl];
	
	for(n=0;n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	for(r=0;r<numk;++r)
	field[n][m][q][r]=new int[numh];
	
	for(n=0;n<numi;++n)
	for(m=0;m<numj;++m)
	for(q=0;q<numk;++q)
	for(r=0;r<numl;++r)
	for(s=0;s<numh;++s)
	field[n][m][q][r][s]=0;
}

void resize_class::Iarray(int **& field, int numi, int *numj)
{
	int n,m;
	
	if(numi>0)
	field=new int*[numi];

	for(n=0;n<numi;++n)
	field[n]=new int[numj[n]];

	for(n=0; n<numi;++n)
	for(m=0; m<numj[n];++m)
	field[n][m]=0;
}

void resize_class::Iarray(int ***& field, int numi, int *numj, int numk)
{
	int n,m,q;
	
	if(numi>0)
	field=new int**[numi];
	
	for(n=0; n<numi; ++n)
	field[n] = new int*[numj[n]];

	for(n=0; n<numi; ++n)
	for(m=0; m<numj[n]; ++m)
	field[n][m] = new int[numk];

	for(n=0; n<numi;++n)
	for(m=0;m<numj[n];++m)
	for(q=0;q<numk;++q)
	field[n][m][q]=0;
}

void resize_class::cvec(char *& field, int numi)
{	
	if(numi>0)
	field = new char[numi];
}

void resize_class::cvec(char **& field, int numi, int numj)
{
	int i;
	
	if(numi>0)
	field=new char*[numi];

	for(i=0;i<numi;++i)
	field[i]=new char[numj];
}

void resize_class::cvec(char ***& field, int numi, int numj, int numk)
{

	int i,j;
	
	if(numi>0)
	field = new char**[numi];

	for(i=0; i<numi; ++i)
	field[i] = new char*[numj];

	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	field[i][j] = new char[numk];
}

void resize_class::rank(int r)
{
	
	pararank=r;
}
