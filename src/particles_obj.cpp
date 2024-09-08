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
                entries(tracers_obj::entries+(individuals?4+3+1+3+3+3+3:0)), // update when adding more data
                d50(_d50), density(_density)
                
{	
    if(_capacity>0 && entries>tracers_obj::entries)
    {
        U = new double[_capacity];
        V = new double[_capacity];
        W = new double[_capacity];

        ParcelFactor = new double[_capacity];
        
        XRK1 = new double[_capacity];
        YRK1 = new double[_capacity];
        ZRK1 = new double[_capacity];
        
        URK1 = new double[_capacity];
        VRK1 = new double[_capacity];
        WRK1 = new double[_capacity];
        
        Uf = new double[_capacity];
        Vf = new double[_capacity];
        Wf = new double[_capacity];
        shear_eff = new double[_capacity];
        shear_crit = new double[_capacity];
        drag = new double[_capacity];

        fill_data(0,_size);
    }
    else
    {
        U=nullptr;
        V=nullptr;
        W=nullptr;
        ParcelFactor=nullptr;
    }
}

/// @brief Deletes individual data if existent
particles_obj::~particles_obj()
{
    if(entries>tracers_obj::entries)
    {
        delete[] U;
        U=nullptr;
        delete[] V;
        V=nullptr;
        delete[] W;
        W=nullptr;

        delete[] XRK1;
        XRK1=nullptr;
        delete[] YRK1;
        YRK1=nullptr;
        delete[] ZRK1;
        ZRK1=nullptr;
        
        delete[] URK1;
        URK1=nullptr;
        delete[] VRK1;
        VRK1=nullptr;
        delete[] WRK1;
        WRK1=nullptr;

        delete[] ParcelFactor;
        ParcelFactor=nullptr;

        delete[] Uf;
        Uf=nullptr;
        delete[] Vf;
        Vf=nullptr;
        delete[] Wf;
        Wf=nullptr;
        delete[] shear_eff;
        shear_eff=nullptr;
        delete[] shear_crit;
        shear_crit=nullptr;
        delete[] drag;
        drag=nullptr;
    }
}

/// \copydoc tracers_obj::debug
void particles_obj::debug()
{
    std::cout<<"particle_obj::debug"<<std::endl;
    // insert code for debugging here //
    
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
        
        double* newXRK1=new double[capacity_desired];
        std::memcpy( newXRK1, XRK1, size * sizeof(double) );
        delete[] XRK1;
        XRK1=newXRK1;

        double* newYRK1=new double[capacity_desired];
        std::memcpy( newYRK1, YRK1, size * sizeof(double) );
        delete[] YRK1;
        YRK1=newYRK1;

        double* newZRK1=new double[capacity_desired];
        std::memcpy( newZRK1, ZRK1, size * sizeof(double) );
        delete[] ZRK1;
        ZRK1=newZRK1;
        
        
        double* newURK1=new double[capacity_desired];
        std::memcpy( newURK1, URK1, size * sizeof(double) );
        delete[] URK1;
        URK1=newURK1;

        double* newVRK1=new double[capacity_desired];
        std::memcpy( newVRK1, VRK1, size * sizeof(double) );
        delete[] VRK1;
        VRK1=newVRK1;

        double* newWRK1=new double[capacity_desired];
        std::memcpy( newWRK1, WRK1, size * sizeof(double) );
        delete[] WRK1;
        WRK1=newWRK1;


        double* newParcelFactor=new double[capacity_desired];
        std::memcpy( newParcelFactor, ParcelFactor, size * sizeof(double) );
        delete[] ParcelFactor;
        ParcelFactor=newParcelFactor;

        double* newUf=new double[capacity_desired];
        std::memcpy( newUf, Uf, size * sizeof(double) );
        delete[] Uf;
        Uf=newUf;

        double* newVf=new double[capacity_desired];
        std::memcpy( newVf, Vf, size * sizeof(double) );
        delete[] Vf;
        Vf=newVf;

        double* newWf=new double[capacity_desired];
        std::memcpy( newWf, Wf, size * sizeof(double) );
        delete[] Wf;
        Wf=newWf;

        double* newshear_eff=new double[capacity_desired];
        std::memcpy( newshear_eff, shear_eff, size * sizeof(double) );
        delete[] shear_eff;
        shear_eff=newshear_eff;

        double* newshear_crit=new double[capacity_desired];
        std::memcpy( newshear_crit, shear_crit, size * sizeof(double) );
        delete[] shear_crit;
        shear_crit=newshear_crit;

        double* newdrag=new double[capacity_desired];
        std::memcpy( newdrag, drag, size * sizeof(double) );
        delete[] drag;
        drag=newdrag;

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

        std::memmove(&ParcelFactor[des],&ParcelFactor[src],sizeof(double)*len);
    }
}

/// @brief Set default data for velocities and packing factor
/// @param start First element to set
/// @param end Element after the last set element
void particles_obj::fill_data(size_t start, size_t end)
{
    for(size_t n=start;n<end;n++)
    {
        U[n]=0.0;
        V[n]=0.0;
        W[n]=0.0;
        
        XRK1[n]=0.0;
        YRK1[n]=0.0;
        ZRK1[n]=0.0;
        
        URK1[n]=0.0;
        VRK1[n]=0.0;
        WRK1[n]=0.0;
        
        ParcelFactor[n]=1;

        Uf[n]=0.0;
        Vf[n]=0.0;
        Wf[n]=0.0;
        shear_eff[n]=0.0;
        shear_crit[n]=0.0;
        drag[n]=0.0;
    }
}