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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"tracers_obj.h"
#include<iostream>

/*
Dangers when using:
size_t overflow when adding something to an object at capacity
*/


tracers_obj::tracers_obj(size_t capacity, size_t size, double scale_factor): scale_factor(scale_factor), flag_no_data(-1), entries(3)
{	
    if(capacity>0)
    {
        if(size>capacity)
            capacity=size;

        X=new double[capacity]; // default value: NULL
        Y=new double[capacity]; // default value: NULL
        Z=new double[capacity]; // default value: NULL

        Flag=new int[capacity]; // default value: -1

        Empty=new size_t[capacity]; // default value: size_t_max-1
        this->capacity=capacity;

        this->size=0;
        this->empty_itr=0;
        this->loopindex=0;
        fill(size);
        
    }
}

tracers_obj::~tracers_obj()
{
    clear();
}



void tracers_obj::debug()
{
    std::cout<<"tracers_obj::debug"<<std::endl;
    // insert code for debugging here //
    
}

void tracers_obj::erase(size_t index)
{
    X[index]=NULL;
    Y[index]=NULL;
    Z[index]=NULL;

    Flag[index]=flag_no_data;

    Empty[++empty_itr]=index;
    --size;
}

void tracers_obj::add(double x, double y, double z, int flag)
{
    X[Empty[empty_itr]]=x;
    Y[Empty[empty_itr]]=y;
    Z[Empty[empty_itr]]=z;

    Flag[Empty[empty_itr]]=flag;

    Empty[empty_itr--]=-1;
    if(!(loopindex>size++))
        loopindex++;
}

void tracers_obj::reserve(size_t capacity_desired)
{
    if(0==capacity_desired)
        capacity_desired=ceil(scale_factor*capacity);
    if (capacity_desired>capacity)
    {
        if (capacity_desired>SIZE_MAX)
            std::__throw_length_error("tracers_obj - max capacity reached");

        double* newX=new double[capacity_desired];
        memcpy( newX, X, size * sizeof(double) );
        delete [] X;
        X=newX;

        double* newY=new double[capacity_desired];
        memcpy( newY, Y, size * sizeof(double) );
        delete [] Y;
        Y=newY;

        double* newZ=new double[capacity_desired];
        memcpy( newZ, Z, size * sizeof(double) );
        delete [] Z;
        Z=newZ;


        int* newFlag=new int[capacity_desired];
        memcpy( newFlag, Flag, size * sizeof(int) );
        delete [] Flag;
        Flag=newFlag;

        size_t* newEmpty=new size_t[capacity_desired];
        memcpy( newEmpty, Empty, size * sizeof(int) );
        delete [] Empty;
        Empty=newEmpty;

        capacity=capacity_desired;
        fill_empty();
    }
}

void tracers_obj::clear()
{
	if(capacity>0)
	{
        delete[] X;
        delete[] Y;
        delete[] Z;

        delete[] Flag;
        delete[] Empty;
    }
}

void tracers_obj::fill(size_t index, bool do_empty)
{
    for(size_t n=size; n<index;++n)
    {
        X[n]=NULL;
        Y[n]=NULL;
        Z[n]=NULL;

        Flag[n]=flag_no_data;
    }
    size=index;
    loopindex=index;
    if(do_empty)
    fill_empty();
}

void tracers_obj::fill_empty()
{
    for(size_t n=empty_itr;n<capacity-size;n++)
    {
        Empty[n]=capacity-n-1;
    }
    if(capacity==size)
        empty_itr=capacity;
    else
        empty_itr=capacity-size-1;
}

bool tracers_obj::check_state(bool first)
{
    if(empty_itr>capacity||-size!=empty_itr||size>capacity)
    {
        if(first)
        {
            fix_state();
            return check_state(false);
        }
        else
            return false;
    }
    else
        return true;
}

void tracers_obj::fix_state()
{
    if(this->size>this->capacity)
    {
        size_t real_size=capacity;
        size_t old_size=size;
        reserve(ceil(this->scale_factor*size));
        this->size=real_size;
        fill(old_size,false);

        this->empty_itr=0;
        size_t temp_Empty[capacity];
        for(size_t n=0;n<real_size;n++)
            if(Empty[n]!=NULL)
                temp_Empty[empty_itr++]=n;
        delete [] Empty;
        this->Empty=temp_Empty;
        fill_empty();
    }
    else
    {
        this->empty_itr=0;
        size_t temp_Empty[capacity];
        delete [] Empty;
        this->Empty=temp_Empty;
        for(size_t n=capacity-1; n>0;--n)
            if(X[n]==NULL&&Y[n]==NULL&&Z[n]==NULL&&Flag[n]==-1)
                Empty[empty_itr++]=n;
    }
}

void tracers_obj::optimize()
{
    // Could be optimized to look ahead if a large section is empty to only do one move operation
    size_t loopchange=0;
    for(int n=empty_itr; n>=0;--n)
        if(Empty[n]<loopindex)
        {
            memorymove(Empty[n],Empty[n]+1,sizeof(double)*(loopindex-Empty[n]+1));
            loopchange++;
        }
    loopindex -= loopchange;
}

void tracers_obj::memorymove(size_t des, size_t src, size_t len)
{
    std::memmove(&X[des],&X[src],len);
    std::memmove(&Y[des],&Y[src],len);
    std::memmove(&Z[des],&Z[src],len);

    std::memmove(&Flag[des],&Flag[src],len);
}

void tracers_obj::add_obj(tracers_obj obj)
{
    if(obj.loopindex>ceil(scale_factor*obj.size))
        obj.optimize();
    if(size+obj.size>capacity)
        reserve(size+obj.size);
    for(size_t n=0;n<obj.loopindex;n++)
        add(obj.X[n],obj.Y[n],obj.Z[n],obj.Flag[n]);
}