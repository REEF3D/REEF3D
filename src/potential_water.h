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

#include"potential.h"
#include"increment.h"
#include"fieldint4.h"

class field;

using namespace std;

#ifndef POTENTIAL_WATER_H_
#define POTENTIAL_WATER_H_

class potential_water : public potential, public increment
{

public:
	potential_water(lexer* p);
	virtual ~potential_water();

	virtual void start(lexer*,fdm*, solver*, ghostcell* pgc);


private:
    void rhs(lexer*,fdm*);
	void ucalc(lexer*,fdm*,field&);
	void vcalc(lexer*,fdm*,field&);
	void wcalc(lexer*,fdm*,field&);
    
    void laplace(lexer*,fdm*,field&);
    void ini_bc(lexer*,fdm*,ghostcell*);
    
    
	double starttime,endtime;
	int count;
	int gcval_pot;
    
    fieldint4 bc;
    
    const double eps;
};

#endif

