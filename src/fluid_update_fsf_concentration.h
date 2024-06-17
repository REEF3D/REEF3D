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

#ifndef FLUID_UPDATE_FSF_CONCENTRATION_H_
#define FLUID_UPDATE_FSF_CONCENTRATION_H_

#include"fluid_update.h"
#include"increment.h"

class fdm;
class lexer;
class ghostcell;
class concentration;

using namespace std;

class fluid_update_fsf_concentration : public fluid_update, increment
{
public:
    fluid_update_fsf_concentration(lexer*, fdm*, ghostcell*, concentration*&);
	virtual ~fluid_update_fsf_concentration();

	virtual void start(lexer*, fdm*, ghostcell*);

private:

    static int iocheck,iter;
    int gcval_ro,gcval_visc;
	const double dx;
    double epsi;
	double visc_air,visc_water,ro_air,ro_water;
	double ro_conc, visc_conc;

	concentration *pconcentration;

};

#endif


