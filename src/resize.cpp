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

void resize_class::Darray(double*& field, int numi)
{
    if(numi > 0)
    {
        field = new double[numi] {0.0};
    }
}

void resize_class::Darray(double**& field, int numi, int numj)
{    
    if(numi > 0)
    {
        field = new double*[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new double[numj] {0.0};
        }
    }
}

void resize_class::Darray(double***& field, int numi, int numj, int numk)
{
    if(numi > 0)
    {
        field = new double**[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new double*[numj];
            
            for(int m = 0; m < numj; ++m)
            {
                field[n][m] = new double[numk] {0.0};
            }
        }
    }
}

void resize_class::Darray(double****& field, int numi, int numj, int numk, int numl)
{
    if(numi>0)
    {
        field = new double***[numi];

        for(int n = 0; n<numi; ++n)
        {
            field[n] = new double**[numj];

            for(int m = 0; m < numj; ++m)
            {
                field[n][m] = new double*[numk];
            
                for(int q = 0; q < numk; ++q)
                {
                    field[n][m][q] = new double[numl] {0.0};
                }
            }
        }
    }
}

void resize_class::Darray(double**& field, int numi, int* numj)
{
    if(numi > 0)
    {
        field = new double*[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new double[numj[n]] {0.0};
        }
    }
}

void resize_class::Iarray(int*& field, int numi)
{
    if(numi>0)
    {
        field = new int[numi] {0};
    }
}

void resize_class::Iarray(int**& field, int numi, int numj)
{
    if(numi>0)
    {
        field = new int*[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new int[numj] {0};
        }
    }
}

void resize_class::Iarray(int***& field, int numi, int numj, int numk)
{
    if(numi > 0)
    {
        field = new int**[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new int*[numj];
            
            for(int m = 0; m < numj; ++m)
            {
                field[n][m] = new int[numk] {0};
            }
        }
    }
}

void resize_class::Iarray(int****& field, int numi, int numj, int numk, int numl)
{
    if(numi > 0)
    {
        field = new int***[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new int**[numj];

            for(int m = 0; m < numj; ++m)
            {
                field[n][m] = new int*[numk];
            
                for(int q = 0; q < numk; ++q)
                {
                    field[n][m][q] = new int[numl] {0};
                }
            }
        }
    }
}

void resize_class::Iarray(int*****& field, int numi, int numj, int numk, int numl, int numh)
{
    if(numi > 0)
    {
        field = new int****[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new int***[numj];

            for(int m = 0; m < numj; ++m)
            {
                field[n][m] = new int**[numk];
            
                for(int q = 0; q < numk; ++q)
                {
                    field[n][m][q] = new int*[numl];
                    
                    for(int r = 0; r < numl; ++r)
                    {
                        field[n][m][q][r] = new int[numh] {0};
                    }
                }
            }
        }
    }
}

void resize_class::Iarray(int**& field, int numi, int* numj)
{
    if(numi > 0)
    {
        field = new int*[numi];

        for(int n = 0; n < numi; ++n)
        {
            field[n] = new int[numj[n]] {0};
        }
    }
}

void resize_class::Iarray(int***& field, int numi, int* numj, int numk)
{    
    if(numi > 0)
    {
        field = new int**[numi];
        
        for(int n = 0; n < numi; ++n)
        {
            field[n] = new int*[numj[n]];

            for(int m = 0; m < numj[n]; ++m)
            {
                field[n][m] = new int[numk] {0};
            }
        }
    }
}
