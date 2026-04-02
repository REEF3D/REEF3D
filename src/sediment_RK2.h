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

#ifndef SEDIMENT_RK2_H_
#define SEDIMENT_RK2_H_

#include"sediment_f.h"

using namespace std;

class sediment_RK2 : public sediment_f, public increment
{
public:
    sediment_RK2(lexer*,ghostcell*,turbulence*, patchBC_interface*);
	virtual ~sediment_RK2();
    
    // NHFLOW interface

    void RK2_step1_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override {};
    void RK2_step2_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override {};

private:
    sediment_RK2 *psed;
   
	
};

#endif
