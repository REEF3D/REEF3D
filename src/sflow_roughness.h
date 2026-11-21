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

#ifndef SFLOW_ROUGHNESS_H_
#define SFLOW_ROUGHNESS_H_

class lexer;
class fdm2D;
class slice;

using namespace std;

class sflow_roughness 
{

public:

    virtual ~sflow_roughness() = default;

    // bed friction of a stage of the explicit integrators, implicit in the momentum: the factor 1/(1 + dt gamma) of
    // the Euler step of the stage, gamma = c_f |u|/h of the stage state
    // args: factor (cell centres), U, V, WL (stage state)
    virtual void factor(lexer*, fdm2D*, slice&, slice&, slice&, slice&)=0;

};

#endif
