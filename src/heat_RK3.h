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

#ifndef HEAT_RK3_H_
#define HEAT_RK3_H_

#include"heat_print.h"
#include"field4.h"
#include"bcheat.h"

class heat;
using namespace std;

class heat_RK3 :public bcheat, public heat_print
{
public:
    heat_RK3(lexer *, fdm*, ghostcell*,heat*&);
	virtual ~heat_RK3();
    
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);
    virtual void diff_update(lexer*, fdm*, ghostcell*);
    
    field4 thermdiff;
    field4 ark1,ark2,Tdiff;

private:
    void clearrhs(lexer*,fdm*,ghostcell*);
    
	int gcval_heat;
	double starttime, endtime;
};

#endif
