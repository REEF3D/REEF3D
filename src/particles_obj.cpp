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

#include"particles_obj.h"
#include"lexer.h"
#include <cstdint>
#include <cstring>
#include<iostream>

/*
Dangers when using:
size_t overflow when adding something to an object at capacity
*/

/// @brief Container for particles with density, d50 and x,y,z velocities
/// @param _capacity Desired initial capacity
/// @param _d50 \copydoc particles_obj::d50
/// @param _density \copydoc particles_obj::density
/// @param individuals Whether particles are allowed to have individual data besides position
/// @param _size Desired number of partices at default position (0,0,0|INT32_MIN) with individual data (0,0,0|0)
/// @param _scale_factor Sets ::scale_factor for ::reserve
particles_obj::particles_obj(size_t _capacity, double _d50, double _density, bool individuals, size_t _size, double _scale_factor):
                tracers_obj(_capacity,_size,_scale_factor),
                entries(tracers_obj::entries+(individuals?4:0)), // update when adding more data
                d50(_d50), density(_density)
                
{	
    if(_capacity>0 && entries>tracers_obj::entries)
    {
        U = new double[_capacity];
        V = new double[_capacity];
        W = new double[_capacity];
        
        PackingFactor = new double[_capacity];

        fill_data(0,_size);
    }
}

/// @brief Deletes individual data if existent
particles_obj::~particles_obj()
{
    if(entries>tracers_obj::entries)
    {
        delete[] U;
        U=0;
        delete[] V;
        V=0;
        delete[] W;
        W=0;

        delete[] PackingFactor;
        PackingFactor=0;
    }
}

/// \copydoc tracers_obj::debug
void particles_obj::debug()
{
    std::cout<<"particle_obj::debug"<<std::endl;
    // insert code for debugging here //
    
}

/// \copydoc tracers_obj::erase
void particles_obj::erase(size_t index)
{
    tracers_obj::erase(index);

    if(entries>tracers_obj::entries)
    {
        U[index]=0;
        V[index]=0;
        W[index]=0;

        PackingFactor[index]=0;
    }
}

/// \copydoc tracers_obj::erase_all
void particles_obj::erase_all()
{
    tracers_obj::erase_all();

    if(entries>tracers_obj::entries)
    {
        delete[] U;
        delete[] V;
        delete[] W;
        U = new double[capacity];
        V = new double[capacity];
        W = new double[capacity];

        delete[] PackingFactor;
        PackingFactor = new double[capacity];
    }
}

/// @brief Addes new particle with prescribed position and state
/// @param x Position in x-dir
/// @param y Position in y-dir
/// @param z Position in z-dir
/// @param flag State - stationary, moving, etc.
/// @param u Velocity in x-dir
/// @param v Velocity in y-dir
/// @param w Velocity in z-dir
/// @param packingFactor Number of real particles represented by the element
/// @return Index of added particle
size_t particles_obj::add(double x, double y, double z, int flag, double u, double v, double w, double packingFactor)
{
    size_t index=tracers_obj::add(x,y,z,flag);
    if(entries>tracers_obj::entries)
        add_data(index,u,v,w,packingFactor);
    return index;
}

/// \copydoc tracers_obj::reserve
size_t particles_obj::reserve(size_t capacity_desired)
{
    if(0==capacity_desired)
        capacity_desired=ceil(scale_factor*capacity);
    if (capacity_desired>capacity)
    {
        if (capacity_desired>SIZE_MAX)
            std::__throw_length_error("particles_obj - max capacity reached");

        tracers_obj::reserve(capacity_desired);

        double* newU=new double[capacity_desired];
        std::memcpy( newU, U, size * sizeof(double) );
        delete[] U;
        U=newU;

        double* newV=new double[capacity_desired];
        std::memcpy( newV, V, size * sizeof(double) );
        delete[] V;
        V=newV;

        double* newW=new double[capacity_desired];
        std::memcpy( newW, W, size * sizeof(double) );
        delete[] W;
        W=newW;

        double* newPackingFactor=new double[capacity_desired];
        std::memcpy( newPackingFactor, PackingFactor, size * sizeof(double) );
        delete[] PackingFactor;
        PackingFactor=newPackingFactor;
    }
    return capacity;
}

/// \copydoc tracers_obj::fill
void particles_obj::fill(size_t index, bool do_empty, int flag)
{
    if(entries>tracers_obj::entries)
        fill_data(size,index);
    tracers_obj::fill(index,do_empty,flag);
}

/// \copydoc tracers_obj::check_state
bool particles_obj::check_state(bool first)
{
    if(capacity-size!=empty_itr||size>capacity)
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

/// \copydoc tracers_obj::fix_state
void particles_obj::fix_state() // ToDo - update
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

/// \copydoc tracers_obj::optimize
void particles_obj::optimize() // ToDo - update
{
    // Could be optimized to look ahead if a large section is empty to only do one move operation
    if(loopindex>(1.0+overhead)*size)
    {
        size_t loopchange=0;
        for(int n=empty_itr; n>=0;--n)
            if(Empty[n]<loopindex)
            {
                memorymove(Empty[n],Empty[n]+1,(loopindex-Empty[n]+1));

                loopchange++;
            }
        loopindex -= loopchange;
        if(loopindex>0)
        reset_Empty();
    }
}

/// \copydoc tracers::memorymove
void particles_obj::memorymove(size_t des, size_t src, size_t len)
{
    tracers_obj::memorymove(des,src,len);

    if(entries>tracers_obj::entries)
    {
        std::memmove(&U[des],&U[src],sizeof(double)*len);
        std::memmove(&V[des],&V[src],sizeof(double)*len);
        std::memmove(&W[des],&W[src],sizeof(double)*len);

        std::memmove(&PackingFactor[des],&PackingFactor[src],sizeof(double)*len);
    }
}

/// \copydoc tracers::add_obj
void particles_obj::add_obj(particles_obj* obj)
{
    if(obj->size>0)
    {
        if(size+obj->size>capacity)
            reserve(size+obj->size);
        
        if(obj->entries>obj->tracers_obj::entries && entries>tracers_obj::entries)
            for(size_t n=0;n<obj->loopindex;n++)
                add(obj->X[n],obj->Y[n],obj->Z[n],obj->Flag[n],obj->U[n],obj->V[n],obj->W[n],obj->PackingFactor[n]);
        else
            tracers_obj::add_obj(obj);
    }
}

/// @brief Additional data input
/// @param index Index of accessed element
/// @param u Velocity in x-dir
/// @param v Velocity in y-dir
/// @param w Velocity in z-dir
/// @param packingFactor Number of real particles represented by the element
void particles_obj::add_data(size_t index, double u, double v, double w, double packingFactor)
{
    U[index] = u;
    V[index] = v;
    W[index] = w;
    PackingFactor[index] = packingFactor;
}

/// @brief Set default data for velocities and packing factor
/// @param start First element to set
/// @param end Element after the last set element
void particles_obj::fill_data(size_t start, size_t end)
{
    for(size_t n=start;n<end;n++)
    {
        U[n]=0;
        V[n]=0;
        W[n]=0;
        PackingFactor[n]=0;
    }
}