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

#include "tracers_obj.h"

#include <math.h>
#include <cstring>
#include <iostream>

/// @brief Tracer style particles\n
/// Contains and manages massless particles used as tracers in the flow field
/// @param capacity Desired initial capacity
/// @param size Desired number of tracers at default position (0,0,0|INT32_MIN)
/// @param scale_factor Sets ::scale_factor for ::reserve
tracers_obj::tracers_obj(size_t _capacity, size_t _size, double _scale_factor): scale_factor(_scale_factor), entries(3)
{	
    if(_capacity>0)
    {
        if(_size>_capacity)
            _capacity=_size;

        X = new double[_capacity]; // default value: NULL
        Y = new double[_capacity]; // default value: NULL
        Z = new double[_capacity]; // default value: NULL

        Flag = new int[_capacity]; // default value: INT32_MIN

        Empty = new size_t[_capacity];
        capacity=_capacity;

        size=0;
        empty_itr=0;
        loopindex=0;
        fill(_size);        
    }
}

/// @brief Deletes the constructed arrays.
tracers_obj::~tracers_obj()
{
    delete[] X;
    X=0;
    delete[] Y;
    Y=0;
    delete[] Z;
    Z=0;

    delete[] Flag;
    Flag=0;
    
    delete[] Empty;
    Empty=0;
}

/// @brief Contains debugging code
void tracers_obj::debug()
{
    std::cout<<"tracers_obj::debug"<<std::endl;
    // insert code for debugging here //
    
}

/// @brief Removes \p _index
void tracers_obj::erase(size_t _index)
{
    X[_index]=NULL;
    Y[_index]=NULL;
    Z[_index]=NULL;

    Flag[_index]=INT32_MIN;

    Empty[++empty_itr]=_index;
    --size;
}

/// @brief Clears all data
void tracers_obj::erase_all()
{
    delete[] X;
    delete[] Y;
    delete[] Z;
    X=new double[capacity];
    Y=new double[capacity];
    Z=new double[capacity];

    delete[] Flag;
    Flag=new int[capacity];

    delete[] Empty;
    Empty=new size_t[capacity];

    size=0;
    empty_itr=0;
    Empty[empty_itr]=0;
    loopindex=0;
    fill_empty();
}

/// @brief Addes new particle with prescribed position and state
/// @param x Position in x-dir
/// @param y Position in y-dir
/// @param z Position in z-dir
/// @param flag State - stationary, moving, etc.
/// @return Index of added particle
size_t tracers_obj::add(double x, double y, double z, int flag)
{
    size_t index=Empty[empty_itr];
    X[index]=x;
    Y[index]=y;
    Z[index]=z;

    Flag[index]=flag;

    Empty[empty_itr--]=-1;
    if(!(loopindex>size++))
        loopindex++;
    return index;
}

/// @brief Addes new element based on \p _index of \p obj
/// @return Index of newly added element
size_t tracers_obj::add_entry(tracers_obj* obj, size_t _index)
{
    size_t index=Empty[empty_itr];
    X[index]=obj->X[_index];
    Y[index]=obj->Y[_index];
    Z[index]=obj->Z[_index];

    Flag[index]=obj->Flag[_index];

    Empty[empty_itr--]=-1;
    if(!(loopindex>size++))
        loopindex++;
    return index;
}

/// @brief Reserves memory for new capacity \p capacity_desired
/// @return Actual new capacity
size_t tracers_obj::reserve(size_t capacity_desired)
{
    if(0==capacity_desired)
        capacity_desired=ceil(scale_factor*capacity);
    capacity_desired=ceil(scale_factor*capacity_desired);
    if (capacity_desired>capacity)
    {
        if (capacity_desired>SIZE_MAX)
            std::__throw_length_error("tracers_obj - max capacity reached");

        double* newX=new double[capacity_desired];
        std::memcpy( newX, X, loopindex * sizeof(double) );
        delete[] X;
        X=newX;

        double* newY=new double[capacity_desired];
        std::memcpy( newY, Y, loopindex * sizeof(double) );
        delete[] Y;
        Y=newY;

        double* newZ=new double[capacity_desired];
        std::memcpy( newZ, Z, loopindex * sizeof(double) );
        delete[] Z;
        Z=newZ;


        int* newFlag=new int[capacity_desired];
        std::memcpy( newFlag, Flag, loopindex * sizeof(int) );
        delete[] Flag;
        Flag=newFlag;

        size_t* newEmpty=new size_t[capacity_desired];
        std::memcpy( newEmpty, Empty, loopindex * sizeof(size_t) );
        delete[] Empty;
        Empty=newEmpty;
        
        capacity=capacity_desired;
        fill_empty();
    }
    return capacity;
}

/// @brief Fills with default
/// @param _index to which to fill
/// @param do_empty toggle for ::fill_empty
/// @param flag provide if diffent from default
void tracers_obj::fill(size_t _index, bool do_empty, int flag)
{
    for(size_t n=size; n<_index;++n)
    {
        X[n]=NULL;
        Y[n]=NULL;
        Z[n]=NULL;

        Flag[n]=flag;
    }
    size=_index;
    loopindex=_index;
    if(do_empty)
        fill_empty();
}

/// @brief Fills ::Empty with empty spaces\n
/// Inserts in reverse order into ::Empty
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

/// @brief Checks ::size vs ::capacity vs ::empty_itr
/// @return Safe state
bool tracers_obj::check_state(bool first)
{
    if(empty_itr>capacity||capacity-size!=empty_itr||size>capacity)
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

/// @brief Truncate to capcity if size is over capacity or fix Empty entires
void tracers_obj::fix_state()
{
    if(size>capacity)
    {
        size_t real_size=capacity;
        size_t old_size=size;
        reserve(ceil(scale_factor*size));
        size=real_size;
        fill(old_size,false);

        empty_itr=0;
        size_t temp_Empty[capacity];
        for(size_t n=0;n<real_size;n++)
            if(Empty[n]!=NULL)
                temp_Empty[empty_itr++]=n;
        delete [] Empty;
        Empty=temp_Empty;
        fill_empty();
    }
    else
    {
        empty_itr=0;
        size_t temp_Empty[capacity];
        delete [] Empty;
        Empty=temp_Empty;
        for(size_t n=capacity-1; n>0;--n)
            if(X[n]==NULL&&Y[n]==NULL&&Z[n]==NULL&&Flag[n]==-1)
                Empty[empty_itr++]=n;
    }
}

/// @brief Removes intermediate empties
void tracers_obj::optimize()
{
    // Could be optimized to look ahead if a large section is empty to only do one move operation
    size_t loopchange=0;
    for(int n=empty_itr; n>=0;--n)
        if(Empty[n]<loopindex)
        {
            memorymove(Empty[n],Empty[n]+1,loopindex-Empty[n]+1);
            loopchange++;
        }
    loopindex -= loopchange;
    if(loopchange>0)
    reset_Empty();
}

/// @brief Moves section of data memory
/// @param des Index of desitnation
/// @param src Index of source
/// @param len Number elements to move
void tracers_obj::memorymove(size_t des, size_t src, size_t len)
{
    std::memmove(&X[des],&X[src],sizeof(double)*len);
    std::memmove(&Y[des],&Y[src],sizeof(double)*len);
    std::memmove(&Z[des],&Z[src],sizeof(double)*len);

    std::memmove(&Flag[des],&Flag[src],sizeof(int)*len);
}

/// @brief Adds contens of object of type tracers_obj
/// @param obj Inputed data
void tracers_obj::add_obj(tracers_obj* obj)
{
    if(obj->size>0)
    {
        if(obj->loopindex>ceil(scale_factor*obj->size))
            obj->optimize();
        if(size+obj->size>capacity)
            reserve(size+obj->size);
        for(size_t n=0;n<obj->loopindex;n++)
            add(obj->X[n],obj->Y[n],obj->Z[n],obj->Flag[n]);
    }
}

/// @brief Prints position of element to cout
/// @param _index Element to print
void tracers_obj::print(size_t _index)
{
    std::cout<<"Tracer_obj["<<_index<<"]=("<<X[_index]<<","<<Y[_index]<<","<<Z[_index]<<")"<<std::endl;
}

/// @brief Recreates Empty array
void tracers_obj::reset_Empty()
{
    delete[] Empty;
    Empty=new size_t[capacity];
    empty_itr=0;
    fill_empty();
}

