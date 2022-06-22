/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"suspended.h"
#include"increment.h"#include"field3.h"

class turbulence;

using namespace std;

#ifndef SUSPENDED_RK3_H_
#define SUSPENDED_RK3_H_

class suspended_RK3 : public suspended, public increment
{
public:
	suspended_RK3(lexer *, fdm*,turbulence*);
	virtual ~suspended_RK3();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, sediment_fdm*);
	virtual void ctimesave(lexer*, fdm*);    void suspsource(lexer*,fdm*,field&);    void bcsusp_start(lexer*,fdm*,ghostcell*,field&);	void sedfsf(lexer*,fdm*,field&);	void clearrhs(lexer*,fdm*);

	int gcval_susp;

private:
    double starttime;    void fill_wvel(lexer*,fdm*,ghostcell*,sediment*);     field3 wvel;    int count,q;

};

#endif

