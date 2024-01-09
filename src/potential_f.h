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

#include"potential.h"
#include"increment.h"
#include"fieldint4.h"

class field;

using namespace std;

#ifndef POTENTIAL_F_H_
#define POTENTIAL_F_H_

class potential_f : public potential, public increment
{

public:
	potential_f(lexer*);
	virtual ~potential_f();

	virtual void start(lexer*,fdm*, solver*, ghostcell*);


private:
    void rhs(lexer*,fdm*);
	void ucalc(lexer*,fdm*,field&);
	void vcalc(lexer*,fdm*,field&);
	void wcalc(lexer*,fdm*,field&);
    
    void smoothen(lexer*,fdm*,ghostcell*);
    
    void laplace(lexer*,fdm*,ghostcell*,field&);
    void ini_bc(lexer*,fdm*,ghostcell*);
    
    
	double starttime,endtime;
	int count;
	int gcval_pot;
    
    fieldint4 bc;
};

#endif

