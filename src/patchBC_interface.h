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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"increment.h"

class lexer;
class fdm;
class fdm2D;
class ghostcell;
class field;
class slice;

#include"patch_obj.h"

using namespace std;

#ifndef PATCHBC_INTERFACE_H_
#define PATCHBC_INTERFACE_H_

class patchBC_interface
{
public:
    
    virtual void patchBC_ini(lexer*, ghostcell*)=0;
    
    // BC update
    virtual void patchBC_ioflow(lexer*, fdm*, ghostcell*, field&,field&,field&)=0;
    virtual void patchBC_rkioflow(lexer*, fdm*, ghostcell*, field&,field&,field&)=0;
    virtual void patchBC_discharge(lexer*, fdm*, ghostcell*)=0;
    virtual void patchBC_pressure(lexer*, fdm*, ghostcell*, field&)=0;
    virtual void patchBC_waterlevel(lexer*, fdm*, ghostcell*, field&)=0;
    
    
    virtual void patchBC_ioflow2D(lexer*, ghostcell*, slice&, slice&, slice&, slice&)=0;
    virtual void patchBC_rkioflow2D(lexer*, ghostcell*, slice&, slice&, slice&, slice&)=0;
    virtual void patchBC_discharge2D(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&)=0;
    virtual void patchBC_pressure2D(lexer*, ghostcell*, slice&)=0;
    virtual void patchBC_pressure2D_ugrad(lexer*, fdm2D*, slice&, slice&)=0;
    virtual void patchBC_pressure2D_vgrad(lexer*, fdm2D*, slice&, slice&)=0;
    virtual void patchBC_waterlevel2D(lexer*, fdm2D*, ghostcell*, slice&)=0;
    
    virtual void patchBC_loop2D(lexer*, fdm2D*, int&, int&, int&, int&)=0;
    
    patch_obj **patch;
    int obj_count;
    
};

#endif