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

#include"heat_print.h"
#include"field4.h"
#include"bcheat.h"

class heat;

using namespace std;

#ifndef HEAT_RK2_H_
#define HEAT_RK2_H_

class heat_RK2 :public bcheat, public heat_print
{
public:
    heat_RK2(lexer *, fdm*, ghostcell*,heat*&);
	virtual ~heat_RK2();

	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);
    virtual void diff_update(lexer*, fdm*, ghostcell*);
    
    field4 thermdiff;

private:
    void clearrhs(lexer*,fdm*,ghostcell*);
    
	int gcval_heat;
	double starttime, endtime;
};

#endif
