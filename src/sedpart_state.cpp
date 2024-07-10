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

#include "sedpart.h"
#include "sedpart_movement.h"

/// @brief Write out particle data to state file
/// @param result statefile
void sedpart::write_state_particles(lexer *p, ofstream &result)
{
    float ffn=num;
    result.write((char*)&ffn, sizeof (float));
    ffn=volume0;
    result.write((char*)&ffn, sizeof (float));
    size_t ffs=PP.capacity;
    result.write((char*)&ffs, sizeof (size_t));
    ffs=PP.size;
    result.write((char*)&ffs, sizeof (size_t));
    PARTICLELOOP
    {
        ffn=PP.X[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.Y[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.Z[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.Flag[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.U[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.V[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.W[n];
        result.write((char*)&ffn, sizeof (float));
        ffn=PP.PackingFactor[n];
        result.write((char*)&ffn, sizeof (float));
    }
    movement->writeState(p,result);
}

/// @brief Read in particle data from state file
/// @param result statefile
void sedpart::read_state_particles(lexer *p, ifstream& result)
{
    float ffn;
    result.read((char*)&ffn, sizeof (float));
    printcount=ffn;
    result.read((char*)&ffn, sizeof (float));
    volume0=ffn;
    PP.erase_all();
    size_t ffs;
    result.read((char*)&ffs, sizeof (size_t));
    maxparticle=size_t(ffs);
    PP.reserve(maxparticle);
    result.read((char*)&ffs, sizeof (size_t));
    double x,y,z,flag,u,v,w,packing;
    for(size_t n=0; n<ffs;n++)
    {
        result.read((char*)&ffn, sizeof (float));
        x=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        y=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        z=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        flag=int(ffn);
        result.read((char*)&ffn, sizeof (float));
        u=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        v=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        w=double(ffn);
        result.read((char*)&ffn, sizeof (float));
        packing=double(ffn);
        PP.add(x,y,z,flag,u,v,w,packing);
    }
    movement->readState(p,result);
}