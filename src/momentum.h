/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include<vector>

class lexer;
class fdm;
class ghostcell;
class vrans;
class convection;
class diffusion;
class solver;
class ghostcell;
class ioflow;
class turbulence;
class pressure;
class poisson;
class momentum;
class net;
class sixdof;

using namespace std;

#ifndef MOMENTUM_H_
#define MOMENTUM_H_

class momentum
{
public:
	virtual void start(lexer*, fdm*, ghostcell*, vrans*, sixdof*, vector<net*>&)=0;
    virtual void utimesave(lexer*,fdm*, ghostcell*)=0;
    virtual void vtimesave(lexer*,fdm*, ghostcell*)=0;
    virtual void wtimesave(lexer*,fdm*, ghostcell*)=0;

};

#endif
