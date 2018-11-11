/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
class concentration;

using namespace std;

#ifndef FLUID_UPDATE_FSF_ENTRAIN_H_
#define FLUID_UPDATE_FSF_ENTRAIN_H_

class fluid_update_fsf_entrain : public fluid_update, increment
{
public:
    fluid_update_fsf_entrain(lexer*, fdm*, ghostcell*, concentration*&);
	virtual ~fluid_update_fsf_entrain();

	virtual void start(lexer*, fdm*, ghostcell*);
	virtual void start3(lexer*, fdm*, ghostcell*,field&,field&);

private:

    static int iocheck,iter;
    int gcval_ro,gcval_visc;
	const double dx;
    double epsi;
	double visc_air,visc_water,ro_air,ro_water;
	double ro_conc,visc_conc;

	concentration *pconcentration;

};

#endif


