/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class reinitopo;
class discrete;
class ioflow;

using namespace std;

#ifndef SOLID_H_
#define SOLID_H_

class solid : public increment
{
public:
	solid(lexer*, fdm*, ghostcell*);
	virtual ~solid();
	virtual void start(lexer*, fdm*, ghostcell*, ioflow*, discrete*, reinitopo*);

private:

    void solid_topo(lexer*,fdm*,ghostcell*);
    int conv(double);

    int istart, iend, jstart, jend, kstart, kend;
    int qn;


};

#endif


