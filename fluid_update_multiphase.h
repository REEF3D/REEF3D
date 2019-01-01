/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fluid_update.h"
#include"increment.h"

class fdm;
class lexer;
class ghostcell;

using namespace std;

#ifndef FLUID_UPDATE_MULTIPHASE_H_
#define FLUID_UPDATE_MULTIPHASE_H_

class fluid_update_multiphase : public fluid_update, increment
{
public:
    fluid_update_multiphase(lexer*, fdm*, ghostcell*);
	virtual ~fluid_update_multiphase();

	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void start3(lexer*, fdm*, ghostcell*,field&,field&);

private:
    static int iocheck,iter;
    int gcval_ro,gcval_visc;
	int n;
	const double dx,visc3,visc2,visc1,ro1,ro2,ro3;
	double eps12,eps13,eps23;
    double epsi;
};

#endif


