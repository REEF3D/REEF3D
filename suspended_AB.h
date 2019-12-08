/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"suspended.h"
#include"bcsusp.h"
#include"susprhs.h"
#include"field4.h"

class turbulence;

using namespace std;

#ifndef SUSPENDED_AB_H_
#define SUSPENDED_AB_H_

class suspended_AB : public suspended, public bcsusp, public susprhs
{
public:
	suspended_AB(lexer *, fdm*,turbulence*);
	virtual ~suspended_AB();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*);
	virtual void ctimesave(lexer*, fdm*);

	field4 cab;
	int gcval_susp;
private:
    double starttime;

};

#endif
