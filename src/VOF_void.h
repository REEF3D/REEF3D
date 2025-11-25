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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/


#include"freesurface.h"
#include"gradient.h"
#include"norm_vec.h"
#include"reini.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"interpolation.h"

class heat;
class fluid_update;

using namespace std;

#ifndef VOF_VOID_H
#define VOF_VOID_H

class VOF_void : public freesurface, gradient, norm_vec
{
public:
    VOF_void(lexer*, fdm*, ghostcell*,heat*);
    virtual ~VOF_void();
    void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&) override;
    void update(lexer*,fdm*,ghostcell*,field&) override;
    
private:
    fluid_update *pupdate;
};

#endif
