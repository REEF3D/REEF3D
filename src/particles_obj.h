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
class lexer;

class particles_obj : public tracers_obj
{
public:
    particles_obj(size_t=10, double=0.001, double=2650.0, bool=false, size_t=0, double=1.25);
    ~particles_obj();

    void erase(size_t);
    void erase_all();

    size_t reserve(size_t=0);
    void optimize();

    size_t add(double,double,double,int,double=0,double=0,double=0,double=1); // expand when adding additional data
    void add_obj(particles_obj*);
    
    void fill(size_t,bool=true,int=-1);
    void fill_data(size_t,size_t);

    bool check_state(bool=true);
    void debug();
private:
    void fix_state();
    void memorymove(size_t des, size_t src, size_t len);
    void add_data(size_t,double,double,double,double); // expand when adding additional data

public:
    // --- state data ---

    /// \copydoc tracers_obj::entries
    const size_t entries;
    const double overhead = 0.25; 

    // --- particle data ---
    // ---    general    ---

    /// @brief d50 of particle set
    double d50;
    /// @brief Average density of particle set
    const double density;
    
    // ---   individual  ---

    /// @brief Velocity in x-dir
    double* U;
    /// @brief Velocity in y-dir
    double* V;
    /// @brief Velocity in z-dir
    double* W;
    /// @brief Number of real particles represented by the element
    double* PackingFactor;
};

#endif
