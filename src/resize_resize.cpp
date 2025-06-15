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

#include"resize.h"
#include<algorithm>

void resize_class::Dresize(double*& field, int iold, int inew)
{
    double *cache;
    
    Darray(cache,inew);
    
    for(int i = 0; i < std::min(iold,inew); ++i)
    {
        cache[i] = field[i];
    }
    
    del_Darray(field,iold);
    
    field = cache;
}

void resize_class::Dresize(double**& field, int iold, int inew, int jold, int jnew)
{
    double **cache;
    
    Darray(cache,inew,jnew);
    
    for(int i = 0; i < std::min(iold, inew); ++i)
    {
        for(int j = 0; j < std::min(jold,jnew); ++j)
        {
            cache[i][j] = field[i][j];
        }
    }
    
    del_Darray(field,iold,jold);
    
    field = cache;
}

void resize_class::Dresize(double***& field, int iold, int inew, int jold, int jnew, int kold, int knew)
{
    double ***cache;
    
    Darray(cache,inew,jnew,knew);
    
    for(int i = 0; i < std::min(iold, inew); ++i)
    {
        for(int j = 0; j < std::min(jold,jnew); ++j)
        {
            for(int k = 0; k < std::min(kold,knew); ++k)
            {
                cache[i][j][k] = field[i][j][k];
            }
        }
    }
    
    del_Darray(field,iold,jold,kold);
    
    field = cache;
}

void resize_class::Dresize(double****& field, int iold, int inew, int jold, int jnew, int kold, int knew, int lold, int lnew)
{
    double ****cache;
    
    Darray(cache,inew,jnew,knew,lnew);
    
    for(int i = 0; i < std::min(iold, inew); ++i)
    {
        for(int j = 0; j < std::min(jold,jnew); ++j)
        {
            for(int k = 0; k < std::min(kold,knew); ++k)
            {
                for(int l = 0; l < std::min(lold, lnew); ++l)
                {
                    cache[i][j][k][l] = field[i][j][k][l];
                }
            }
        }
    }
    
    del_Darray(field,iold,jold,kold,lold);
    
    field = cache;
}

void resize_class::Iresize(int*& field, int iold, int inew)
{
    int *cache;

    Iarray(cache,inew);
    
    for(int i = 0; i < std::min(iold,inew); ++i)
    {
        cache[i] = field[i];
    }
    
    del_Iarray(field,iold);
    
    field = cache;
}

void resize_class::Iresize(int**& field, int iold, int inew, int jold, int jnew)
{
    int **cache;
    
    Iarray(cache,inew,jnew);
    
    for(int i = 0; i < std::min(iold,inew); ++i)
    {
        for(int j = 0; j < std::min(jold,jnew); ++j)
        {
            cache[i][j] = field[i][j];
        }
    }
    
    del_Iarray(field,iold,jold);
    
    field = cache;
}

void resize_class::Iresize(int***& field, int iold, int inew, int jold, int jnew, int kold, int knew)
{
    int ***cache;

    Iarray(cache,inew,jnew,knew);

    for(int i = 0; i < std::min(iold,inew); ++i)
    {
        for(int j = 0; j < std::min(jold,jnew); ++j)
        {
            for(int k = 0; k < std::min(kold,knew); ++k)
            {
                cache[i][j][k] = field[i][j][k];
            }
        }
    }

    del_Iarray(field,iold,jold,kold);

    field = cache;
}

void resize_class::Iresize(int***& field, int iold, int inew, int *jold, int *jnew, int kold, int knew)
{
    int ***cache;
    
    Iarray(cache,inew,jnew,knew);
    
    for(int i = 0; i < std::min(iold,inew); ++i)
    {
        for(int j = 0; j < std::min(jold[i],jnew[i]); ++j)
        {
            for(int k = 0; k < std::min(kold,knew); ++k)
            {
                cache[i][j][k] = field[i][j][k];
            }
        }
    }
    
    del_Iarray(field,iold,jold,kold);
    
    field = cache;
}

void resize_class::Iresize(int****& field, int iold, int inew, int jold, int jnew, int kold, int knew, int lold, int lnew)
{
    int ****cache;
    
    Iarray(cache,inew,jnew,knew,lnew);
    
    for(int i = 0; i < std::min(iold,inew); ++i)
    {
        for(int j = 0; j < std::min(jold,jnew); ++j)
        {
                for(int k = 0; k < std::min(kold,knew); ++k)
            {
                for(int l = 0; l < std::min(lold,lnew); ++l)
                {
                    cache[i][j][k][l] = field[i][j][k][l];
                }
            }
        }
    }
    
    del_Iarray(field,iold,jold,kold,lold);
    
    field = cache;
}
