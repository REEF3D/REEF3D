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
#include<iostream>

/*
Dangers when using:
size_t overflow when adding something to an object at capacity
*/


particles_obj::particles_obj(size_t capacity, double d50, double density, double porosity, size_t size, double scale_factor):
                d50(d50), density(density), porosity(porosity), scale_factor(scale_factor), tracers_obj(capacity,size,scale_factor),
                flag_inactive(0), flag_bed(1), flag_bed_load(2), flag_suspended_load(3),
                entries(tracers_obj::entries+0) // update when adding more data
{	
    if(capacity>0)
    {
        if(size>capacity)
            capacity=size;

        // Add new individual data here:
        
    }
}

particles_obj::~particles_obj()
{
    clear();
}



void particles_obj::debug()
{
    std::cout<<"particle_obj::debug"<<std::endl;
    // insert code for debugging here //
    
}

void particles_obj::erase(size_t index)
{
    tracers_obj::erase(index);
}

void particles_obj::add(double x, double y, double z, int flag)
{
    tracers_obj::add(x,y,z,flag);
    if(entries>tracers_obj::entries)
        add_data();
}

void particles_obj::reserve(size_t capacity_desired)
{
    if(0==capacity_desired)
        capacity_desired=ceil(scale_factor*capacity);
    if (capacity_desired>capacity)
    {
        if (capacity_desired>SIZE_MAX)
            std::__throw_length_error("particles_obj - max capacity reached");

        tracers_obj::reserve(capacity_desired);
    }
}

void particles_obj::clear()
{
	if(capacity>0)
	{

    }
    tracers_obj::clear();
}

void particles_obj::fill(size_t index, bool do_empty)
{
    if(entries>tracers_obj::entries)
        for(size_t n=size; n<index;++n)
        {
            
        }
    tracers_obj::fill(index,do_empty);
}

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

void particles_obj::fix_state()
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

void particles_obj::optimize()
{
    // Could be optimized to look ahead if a large section is empty to only do one move operation
    size_t loopchange=0;
    for(int n=empty_itr; n>=0;--n)
        if(Empty[n]<loopindex)
        {
            tracers_obj::memorymove(Empty[n],Empty[n]+1,sizeof(double)*(loopindex-Empty[n]+1));

            loopchange++;
        }
    loopindex -= loopchange;
}
void particles_obj::add_obj(particles_obj obj)
{
    if(size+obj.size>capacity)
        reserve(size+obj.size);
    tracers_obj::add_obj(obj);
    if(entries>tracers_obj::entries)
        for(size_t n=0;n<obj.loopindex;n++)
            add_data();
}
void particles_obj::add_data()
{

}