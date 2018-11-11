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
#include"air_entrainment_ibc.h"
#include"field4.h"

class concentration;

using namespace std;

#ifndef AIR_ENTRAINMENT_IM1_H_
#define AIR_ENTRAINMENT_IM1_H_

class air_entrainment_IM1 :public concentration_io, public air_entrainment_ibc
{
public:
    air_entrainment_IM1(lexer *, fdm*, ghostcell*);
	virtual ~air_entrainment_IM1();

	virtual void start(fdm*, lexer*, discrete*, diffusion*, turbulence*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);

	field4 Cn;

private:
    void timesource(lexer* p, fdm* a, field& fn);
    void clearrhs(lexer*,fdm*);

	int gcval_concentration;
	double starttime, endtime;
};

#endif
