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

#include"multiphase_fluid_update.h"
#include"increment.h"

using namespace std;

#ifndef MULTIPHASE_FLUID_UPDATE_F_H_
#define MULTIPHASE_FLUID_UPDATE_F_H_

class multiphase_fluid_update_f : public multiphase_fluid_update, increment
{
public:
    multiphase_fluid_update_f(lexer*, fdm*, ghostcell*);
	virtual ~multiphase_fluid_update_f();

	virtual void start(lexer*, fdm*, ghostcell*,field&,field&);

private:
    static int iocheck,iter;
    int gcval_ro,gcval_visc;
	int n;
	const double dx,visc3,visc2,visc1,ro1,ro2,ro3;
	double eps12,eps13,eps23;
    double epsi;
};

#endif


