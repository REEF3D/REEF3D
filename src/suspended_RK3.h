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

#ifndef SUSPENDED_RK3_H_
#define SUSPENDED_RK3_H_

#include"suspended.h"
#include"increment.h"
#include"field3.h"

using namespace std;

class suspended_RK3 : public suspended, public increment
{
public:
	suspended_RK3(lexer *, fdm*);
	virtual ~suspended_RK3();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, sediment_fdm*);
	virtual void ctimesave(lexer*, fdm*);

    void suspsource(lexer*,fdm*,field&,sediment_fdm*);
    void bcsusp_start(lexer*,fdm*,ghostcell*,sediment_fdm*,field&);
	void sedfsf(lexer*,fdm*,field&);
	void clearrhs(lexer*,fdm*);
    void fillconc(lexer*,fdm*,sediment_fdm*);

	int gcval_susp;

private:
    double starttime;
    void fill_wvel(lexer*,fdm*,ghostcell*,sediment_fdm*); 
    field3 wvel;

    int count,q;

};

#endif

