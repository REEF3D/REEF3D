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

#include"concentration_io.h"
#include"ibc_concentration.h"
#include"field4.h"

class concentration;

using namespace std;

#ifndef CONCENTRATION_IM1_H_
#define CONCENTRATION_IM1_H_

class concentration_IM1 :public ibc_concentration, public concentration_io
{
public:
    concentration_IM1(lexer *, fdm*, ghostcell*);
	virtual ~concentration_IM1();

	virtual void start(fdm*, lexer*, discrete*, diffusion*, turbulence*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);

	field4 Cn;

private:
    void timesource(lexer* p, fdm* a, field& fn);
    void clearrhs(lexer*,fdm*,ghostcell*);

	int gcval_concentration;
	double starttime, endtime;
};

#endif
