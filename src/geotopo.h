/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class reinitopo;
class ioflow;
class vrans;

using namespace std;

#ifndef GEOTOPO_H_
#define GEOTOPO_H_

class geotopo : public increment
{
public:
	geotopo(lexer*, fdm*, ghostcell*);
	virtual ~geotopo();
	virtual void start(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, vrans*);

private:
    void dat(lexer*,fdm*,ghostcell*);

    int conv(double);

    int istart, iend, jstart, jend, kstart, kend;
    int qn;


};

#endif


