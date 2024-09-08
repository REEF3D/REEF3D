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

/// @brief Addes new particle with prescribed position and state
/// @param x Position in x-dir
/// @param y Position in y-dir
/// @param z Position in z-dir
/// @param flag State - stationary, moving, etc.
/// @param u Velocity in x-dir
/// @param v Velocity in y-dir
/// @param w Velocity in z-dir
/// @param parcelFactor Number of real particles represented by the element
/// @return Index of added particle
size_t particles_obj::add(double x, double y, double z, int flag, double u, double v, double w, 
                            double parcelFactor, double xrk1, double yrk1, double zrk1, 
                            double urk1, double vrk1, double wrk1, double uF, double vF, double wF, 
                            double shearEff, double shearCrit, double _drag)
{
    size_t index=tracers_obj::add(x,y,z,flag);
    if(entries>tracers_obj::entries)
        add_data(index,x,y,z,u,v,w,parcelFactor,xrk1,yrk1,zrk1,urk1,vrk1,wrk1,uF,vF,wF,shearEff,shearCrit,_drag);
    return index;
    
 
}

/// \copydoc tracers::add_obj
void particles_obj::add_obj(particles_obj* obj)
{
    if(obj->size>0)
    {
        if(obj->density!=density)
            std::cerr<<"particles_obj::add_obj - density mismatch"<<std::endl;
        if(obj->d50!=d50)
            std::cerr<<"particles_obj::add_obj - d50 mismatch"<<std::endl;
        if(size+obj->size>capacity)
            reserve(size+obj->size);
        
        if(obj->entries>obj->tracers_obj::entries && entries>tracers_obj::entries)
            for(size_t n=0;n<obj->loopindex;n++)
                add(obj->X[n],obj->Y[n],obj->Z[n],obj->Flag[n],obj->U[n],obj->V[n],obj->W[n],
                obj->ParcelFactor[n],obj->XRK1[n],obj->YRK1[n],obj->ZRK1[n],
                obj->URK1[n],obj->VRK1[n],obj->WRK1[n],obj->Uf[n],obj->Vf[n],obj->Wf[n],
                obj->shear_eff[n],obj->shear_crit[n],obj->drag[n]);
        else
            tracers_obj::add_obj(obj);
    }
}

/// @brief Additional data input
/// @param index Index of accessed element
/// @param u Velocity in x-dir
/// @param v Velocity in y-dir
/// @param w Velocity in z-dir
/// @param parcelFactor Number of real particles represented by the element
void particles_obj::add_data(size_t index, double x, double y, double z, 
                            double u, double v, double w, double parcelFactor, 
                            double xrk1, double yrk1, double zrk1, double urk1, double vrk1, double wrk1, 
                            double uF, double vF, double wF, double shearEff, double shearCrit, double _drag)
{
    XRK1[index] = x;
    YRK1[index] = y;
    ZRK1[index] = z;
    
    U[index] = u;
    V[index] = v;
    W[index] = w;

    XRK1[index] = xrk1;
    YRK1[index] = yrk1;
    ZRK1[index] = zrk1;
    
    URK1[index] = urk1;
    VRK1[index] = vrk1;
    WRK1[index] = wrk1;
    
    ParcelFactor[index] = parcelFactor;

    Uf[index] = uF;
    Vf[index] = vF;
    Wf[index] = wF;
    shear_eff[index] = shearEff;
    shear_crit[index] = shearCrit;
    drag[index] = _drag;
}
