/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SIXDOF_LOAD_H_
#define SIXDOF_LOAD_H_

class lexer;
class sixdof_rigidbody;
class sixdof_geometry;

//  External load model acting on a 6DOF body (solver independent), e.g. the ship module:
//  evaluated with the body state of every stage, before the rigid-body right-hand side.
//
//  add_load: add the load to F[0..5] = X, Y, Z, K, M, N (inertial frame, moments about the
//            centre of gravity)
//  print:    output once per time step (rank 0 decides inside)

class sixdof_load
{
public:
    virtual ~sixdof_load() {}
    
    virtual void add_load(lexer*, const sixdof_rigidbody&, const sixdof_geometry&, double*)=0;
    virtual void print(lexer*) {}
};

#endif
