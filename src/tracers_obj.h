/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef TRACERSOBJ_H_
#define TRACERSOBJ_H_

#include <stdio.h>
#include <math.h>
#include <cstring>

/*
Philosophy: performance, memory usage, ease of use

All data should be stored as at max a 1D array
Internal data access should be via iterator
No external access to iterator

Out of bounds safe
Thread safe
*/

class tracers_obj
{
public:
    tracers_obj(size_t, size_t=0, double=1.25);
    virtual ~tracers_obj();

    void erase(size_t);
    void reserve(size_t=0);
    void add(double,double,double,int);
    bool check_state(bool=true);
    void optimize();
    void debug();
    void add_obj(tracers_obj);

protected:
    void fill(size_t,bool=true);
    void fill_empty();
    void fix_state();
    void clear();
    void memorymove(size_t, size_t, size_t);

public:
    // state data
    size_t size;
    size_t loopindex;
    size_t capacity;
    const size_t entries;
    const int flag_no_data;

    // tracer data
    double* X;
    double* Y;
    double* Z;

    int* Flag;

protected:
    size_t empty_itr;
    size_t* Empty;
    const double scale_factor;
};

#endif
