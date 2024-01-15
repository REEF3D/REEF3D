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

#ifndef PARTICLESOBJ_H_
#define PARTICLESOBJ_H_

#include "tracers_obj.h"

/*
Philosophy: performance, memory usage, ease of use

All data should be stored as at max a 1D array
Internal data access should be via iterator
No external access to iterator

Out of bounds safe
Thread safe
*/

class particles_obj : public tracers_obj
{
public:
    particles_obj(size_t, double=0.001, double=2650.0, double=0.4, size_t=0, double=1.25);
    ~particles_obj();

    void erase(size_t);
    void reserve(size_t=0);
    void add(double,double,double,int); // expand when adding additional data
    void add_obj(particles_obj);
    bool check_state(bool=true);
    void optimize();
    void debug();

private:
    void fill(size_t,bool=true);
    void fix_state();
    void clear();
    void add_data(); // expand when adding additional data

public:
    // state data
    const size_t entries;
    const int flag_inactive;
    const int flag_bed;
    const int flag_bed_load;
    const int flag_suspended_load;

    //particle data
    //general
    double d50;
    const double density;
    const double porosity;

    //individual


private:
    const double scale_factor;
};

#endif
