/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"timestep_ptf.h"
#include"increment.h"

class turbulence;

using namespace std;

#ifndef FIXTIMESTEP_PTF_H_
#define FIXTIMESTEP_PTF_H_

class fixtimestep_ptf : public timestep_ptf, public increment
{
public:
	fixtimestep_ptf(lexer*);
	virtual ~fixtimestep_ptf();
	virtual void start(fdm_ptf*,lexer*,ghostcell*,turbulence*);
	virtual void ini(fdm_ptf*,lexer*,ghostcell*);


};

#endif

