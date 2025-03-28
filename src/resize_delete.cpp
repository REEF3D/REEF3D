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

void resize_class::del_Darray(double*& field, int numi)
{
    if(numi > 0)
    {
        delete [] field;
    }
}


void resize_class::del_Darray(double**& field, int numi, int numj)
{
    if(numi > 0)
    {
        if(numj > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                delete [] field[i];
            }
        }
        delete [] field;
    }
}

void resize_class::del_Darray(double***& field, int numi, int numj, int numk)
{
    if(numi > 0)
    {
        if(numj > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                if(numk > 0)
                {
                    for(int j = 0; j < numj; ++j)
                    {
                        delete [] field[i][j];
                    }
                }
                delete [] field[i];
            }
        }
        delete [] field;
    }
}

void resize_class::del_Darray(double****& field, int numi, int numj, int numk, int numl)
{
    if(numi > 0)
    {
        if(numj > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                if(numk > 0)
                {
                    for(int j = 0; j < numj; ++j)
                    {
                        if(numl > 0)
                        {
                            for(int k = 0; k < numk; ++k)
                            {
                                delete [] field[i][j][k];
                            }
                        }
                        delete [] field[i][j];
                    }
                }
                delete [] field[i];
            }
        }
        delete [] field;
    }
}


void resize_class::del_Iarray(int*& field, int numi)
{    
    if(numi > 0)
    {
        delete [] field;
    }
}


void resize_class::del_Iarray(int**& field, int numi, int numj)
{
    int i;
    
    if(numj>0)
    for(i=0;i<numi;++i)
    delete [] field[i];
    
    if(numi>0)
    delete [] field;
}


void resize_class::del_Iarray(int***& field, int numi, int numj, int numk)
{
    if(numi > 0)
    {
        if(numj > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                delete [] field[i];
            }
        }
        delete [] field;
    }
}

void resize_class::del_Iarray(int***& field, int numi, int* numj, int numk)
{
    if(numi > 0)
    {
        if(numk > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                for(int j = 0; j < numj[i]; ++j)
                {
                    delete [] field[i][j];
                }
                delete [] field[i];
            }
        }
        delete [] field;
    }
}

void resize_class::del_Iarray(int****& field, int numi, int numj, int numk, int numl)
{
    if(numi > 0)
    {
        if(numj > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                if(numk > 0)
                {
                    for(int j = 0; j < numj; ++j)
                    {
                        delete [] field[i][j];
                    }
                }
                delete [] field[i];
            }
        }
        delete [] field;
    }
}

void resize_class::del_Iarray(int*****& field,int numi, int numj, int numk, int numl, int numh)
{
    if(numi > 0)
    {
        if(numj > 0)
        {
            for(int i = 0; i < numi; ++i)
            {
                if(numk > 0)
                {
                    for(int j = 0; j < numj; ++j)
                    {
                        if(numl > 0)
                        {
                            for(int k = 0; k < numk; ++k)
                            {
                                if(numh > 0)
                                {
                                    for(int l = 0; l < numl; ++l)
                                    {
                                        delete[] field[i][j][k][l];
                                    }
                                }
                                delete[] field[i][j][k];
                            }
                        }
                        delete[] field[i][j];
                    }
                }
                delete[] field[i];
            }
        }
        delete[] field;
    }
}
