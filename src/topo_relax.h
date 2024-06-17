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

#ifndef TOPO_RELAX_H_
#define TOPO_RELAX_H_

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class sediment_fdm;
class slice;

using namespace std;

class topo_relax : public increment
{
public:
    topo_relax(lexer*);
    virtual ~topo_relax();

	virtual void start(lexer*,ghostcell*,sediment_fdm*);
    virtual double rf(lexer*,ghostcell*);

private:
	double distcalc(lexer*, double, double, double);
	double r1(lexer*, double, double);
	
	double *tan_betaS73,*betaS73,*dist_S73;
	double val;
	

};

#endif

