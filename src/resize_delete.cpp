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

#include"resize.h"

void resize_class::del_Darray(double *&field,int numi)
{
	if(numi>0)
	delete [ ] field;
}


void resize_class::del_Darray(double **&field,int numi, int numj)
{
	int i;
	
	if(numj>0)
    for(i=0;i<numi;++i)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}

void resize_class::del_Darray(double ***&field,int numi, int numj, int numk)
{
	int i,j;
	
	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	delete [ ] field[i][j];

	if(numj>0)
	for(i=0; i<numi; ++i)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}

void resize_class::del_Darray(double ****&field,int numi, int numj, int numk, int numl)
{
	int i,j,k;
	
	if(numl>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	for(k=0; k<numk; ++k)
	delete [ ] field[i][j][k];

	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	delete [ ] field[i][j];
	
	if(numj>0)
	for(i=0; i<numi; ++i)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}


void resize_class::del_Iarray(int *&field,int numi)
{	
	if(numi>0)
	delete [ ] field;
}


void resize_class::del_Iarray(int **&field,int numi, int numj)
{
	int i;
	
	if(numj>0)
    for(i=0;i<numi;++i)
	delete [ ] field[i];
	
	if(numi>0)
	delete [ ] field;
}


void resize_class::del_Iarray(int ***&field,int numi, int numj, int numk)
{
	int i,j;

	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	delete [ ] field[i][j];

	if(numj>0)
	for(i=0; i<numi; ++i)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}

void resize_class::del_Iarray(int ***&field,int numi, int *numj, int numk)
{
	int i,j;

	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj[i]; ++j)
	delete [ ] field[i][j];

	
	for(i=0; i<numi; ++i)
    if(numj[i]>0)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}

void resize_class::del_Iarray(int ****&field,int numi, int numj, int numk, int numl)
{
	int i,j,k;

	if(numl>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	for(k=0; k<numk; ++k)
	delete [ ] field[i][j][k];

	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	delete [ ] field[i][j];
	
	if(numj>0)
	for(i=0; i<numi; ++i)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}

void resize_class::del_Iarray(int *****&field,int numi, int numj, int numk, int numl, int numh)
{
	int i,j,k,l;

	if(numh>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	for(k=0; k<numk; ++k)
	for(l=0; l<numl; ++l)
	delete [ ] field[i][j][k][l];

	if(numl>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	for(k=0; k<numk; ++k)
	delete [ ] field[i][j][k];
	
	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	delete [ ] field[i][j];
	
	if(numj>0)
	for(i=0; i<numi; ++i)
	delete [ ] field[i];

	if(numi>0)
	delete [ ] field;
}


void resize_class::del_cvec(char *&field,int numi)
{
	if(numi>0)
    delete [ ] field;
}


void resize_class::del_cvec(char **&field,int numi, int numj)
{
	int i;
	
	if(numj>0)
    for(i=0;i<numi;++i)
	delete [ ] field[i];
	
	if(numi>0)
	delete [ ] field;
}

void resize_class::del_cvec(char ***&field,int numi, int numj, int numk)
{
	int i,j;
	
	if(numk>0)
	for(i=0; i<numi; ++i)
	for(j=0; j<numj; ++j)
	delete [ ] field[i][j];
	
	if(numj>0)
	for(i=0; i<numi; ++i)
	delete [ ] field[i];
	
	if(numi>0)
	delete [ ] field;
}

