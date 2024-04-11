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
#include <cstdint>

/*
Philosophy: performance, memory usage, ease of use

All data should be stored as at max a 1D array
Internal data access should be via iterator
No external access to iterator
*/

class tracers_obj
{
public:
    tracers_obj(size_t _capacity=10, size_t _size=0, double _scale_factor=1.25);
    virtual ~tracers_obj();

    virtual void erase(size_t _index);
    void erase_all();
    size_t reserve(size_t capacity_desired=0);
    size_t add(double x, double y, double z, int flag);
    bool check_state(bool first=true);
    void optimize();
    void debug();
    void add_obj(tracers_obj* obj);
    size_t add_entry(tracers_obj* obj, size_t _index);
    void fill(size_t _index, bool do_empty=true, int flag=INT32_MIN);
    void print(size_t _index);

protected:
    void fill_empty();
    void fix_state();
    void memorymove(size_t des, size_t src, size_t len);
    void reset_Empty();

public:
    // --- state data ---

    /// @brief Number of stored particles
    size_t size;
    /// @brief Index for looping\n
    /// ::size is not sufficent if data is sparse
    size_t loopindex;
    /// @brief Currentl capacity to store particles
    size_t capacity;
    /// @brief Number of different data values per particle
    const size_t entries;

    // --- tracer data ---

    /// @brief x-position array
    double* X;
    /// @brief y-position array
    double* Y;
    /// @brief z-position array
    double* Z;

    /// @brief Flags of particles
    int* Flag;

protected:
    /// @brief Current index in ::Empty
    size_t empty_itr;
    /// @brief All empty positions in main data
    size_t* Empty;
    /// @brief Default factor for ::reserve
    const double scale_factor;
};

#endif
