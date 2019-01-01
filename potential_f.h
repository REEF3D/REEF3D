/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"potential.h"
#include"increment.h"

using namespace std;

#ifndef POTENTIAL_F_H_
#define POTENTIAL_F_H_

class potential_f : public potential, public increment
{

public:
	potential_f(lexer* p);
	virtual ~potential_f();

	virtual void start(fdm*,lexer*, solver*, ghostcell* pgc);
	virtual void rhs(lexer*,fdm*);
	virtual void ucalc(lexer*,fdm*);
	virtual void vcalc(lexer*,fdm*);
	virtual void wcalc(lexer*,fdm*);



private:
    void laplace(lexer*,fdm*);
	double starttime,endtime;
	int count;
	int gcval_pot;
};

#endif

