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

#ifndef CONCENTRATION_RK3_H_
#define CONCENTRATION_RK3_H_

#include"concentration_io.h"
#include"field4.h"
#include"bc_concentration.h"

class concentration;
using namespace std;

class concentration_RK3 :public bc_concentration, public concentration_io
{
public:
    concentration_RK3(lexer *, fdm*, ghostcell*);
	virtual ~concentration_RK3();

	virtual void start(fdm*, lexer*, convection*, diffusion*, turbulence*, solver*, ghostcell*, ioflow*);
	virtual void ttimesave(lexer*, fdm*);

private:
    void clearrhs(lexer*,fdm*,ghostcell*);
        
	int gcval_concentration;
	double starttime, endtime;
    
};

#endif
