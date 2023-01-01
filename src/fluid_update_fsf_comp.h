/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fluid_update.h"
#include"increment.h"

class fdm;
class lexer;
class ghostcell;

using namespace std;

#ifndef FLUID_UPDATE_FSF_COMP_H_
#define FLUID_UPDATE_FSF_COMP_H_

class fluid_update_fsf_comp : public fluid_update, increment
{
public:
    fluid_update_fsf_comp(lexer*, fdm*, ghostcell*);
	virtual ~fluid_update_fsf_comp();

	virtual void start(lexer*, fdm*, ghostcell*);

private:
    static int iocheck,iter;
    int gcval_ro,gcval_visc;
	int n;
	double ro_air,epsi;
	const double dx,visc_air,visc_water,ro_water;

};

#endif


