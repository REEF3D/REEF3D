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

class lexer;
class fdm;
class convection;
class ghostcell;
class ioflow;
class reinitopo;
class suspended;
class bedload;
class topo;
class reinitopo;

#include<fstream>

using namespace std;

#ifndef SEDIMENT_INTERFACE_H_
#define SEDIMENT_INTERFACE_H_

class sediment_interface
{
public:

	virtual void start(lexer*, fdm*, convection*, ghostcell*, ioflow*, topo*, reinitopo*, suspended*, bedload*)=0;
	virtual void update(lexer*,fdm*,ghostcell*,ioflow*)=0;
    virtual void relax(lexer*,fdm*,ghostcell*)=0;
    virtual void ini(lexer*,fdm*,ghostcell*)=0;


};

#endif
