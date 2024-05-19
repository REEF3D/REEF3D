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

#include"concentration_io.h"
#include"field4.h"
#include"bc_concentration.h"

class concentration;

using namespace std;

#ifndef CONCENTRATION_AB_H_
#define CONCENTRATION_AB_H_

class concentration_AB :public bc_concentration, public concentration_io
{
public:
    concentration_AB(lexer *, fdm*, ghostcell*);
	virtual ~concentration_AB();

	virtual void start(fdm*, lexer*, convection*, diffusion*, turbulence*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);

private:
    void clearrhs(lexer*,fdm*,ghostcell*);
    
	field4 cab;
	
	int gcval_concentration;
	double starttime, endtime;
};

#endif
