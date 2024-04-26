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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"increment.h"
#include"slice4.h"

class lexer;
class ghostcell;
class fdm_nhf;
class slice;
class patchBC_interface;

#ifndef NHFLOW_RECONSTRUCT_H_
#define NHFLOW_RECONSTRUCT_H_

using namespace std;

class nhflow_reconstruct
{
public:

    virtual void reconstruct_2D_x(lexer*,ghostcell*,fdm_nhf*,slice&,slice&,slice&)=0;
    virtual void reconstruct_2D_y(lexer*,ghostcell*,fdm_nhf*,slice&,slice&,slice&)=0;
    virtual void reconstruct_2D_WL(lexer*,ghostcell*,fdm_nhf*)=0;
    
    virtual void reconstruct_3D_x(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*)=0;
    virtual void reconstruct_3D_y(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*)=0;
    virtual void reconstruct_3D_z(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*)=0;

};

#endif
