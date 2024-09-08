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


/// \copydoc tracers_obj::erase
void particles_obj::erase(size_t index)
{
    tracers_obj::erase(index);

    if(entries>tracers_obj::entries)
    {
        U[index]=0;
        V[index]=0;
        W[index]=0;
        
        URK1[index]=0;
        VRK1[index]=0;
        WRK1[index]=0;

        ParcelFactor[index]=0;

        Uf[index]=0;
        Vf[index]=0;
        Wf[index]=0;
        shear_eff[index]=0;
        shear_crit[index]=0;
        drag[index]=0;
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
        
        delete[] URK1;
        delete[] VRK1;
        delete[] WRK1;
        URK1 = new double[capacity];
        VRK1 = new double[capacity];
        WRK1 = new double[capacity];

        delete[] ParcelFactor;
        ParcelFactor = new double[capacity];

        delete[] Uf;
        delete[] Vf;
        delete[] Wf;
        delete[] shear_eff;
        delete[] shear_crit;
        delete[] drag;
        Uf = new double[capacity];
        Vf = new double[capacity];
        Wf = new double[capacity];
        shear_eff = new double[capacity];
        shear_crit = new double[capacity];
        drag = new double[capacity];
    }
}
