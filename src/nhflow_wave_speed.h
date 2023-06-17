/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#ifndef NHFLOW_WAVE_SPEED_H_
#define NHFLOW_WAVE_SPEED_H_

using namespace std;

class nhflow_wave_speed : public increment
{
public:
	nhflow_wave_speed(lexer*);
	virtual ~nhflow_wave_speed();

    virtual void wave_speed_update(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*);


private:
    nhflow_reconstruct *precon;
};

#endif
