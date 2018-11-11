/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"heat_print.h"
#include"field4.h"
#include"ibcheat.h"
#include"fluid_update.h"

class heat;

using namespace std;

#ifndef HEAT_IM2_H_
#define HEAT_IM2_H_

class heat_IM2 :public ibcheat, public heat_print
{
public:
    heat_IM2(lexer *, fdm*, ghostcell*,heat*&);
	virtual ~heat_IM2();

	virtual void start(fdm*, lexer*, discrete*, diffusion*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);
    virtual void diff_update(lexer*, fdm*, ghostcell*);

	field4 Tn,Tnn;
    field4 thermdiff;

private:
    void timesource(lexer* p, fdm* a, field& fn);
    void clearrhs(lexer*,fdm*,ghostcell*);

	int gcval_heat;
	double starttime, endtime;
	
};

#endif
