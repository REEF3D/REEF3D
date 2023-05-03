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

#include"reinitopo.h"
#include"ddweno.h"
#include"vec.h"
#include"increment.h"

class reinidisc;
class picard;

using namespace std;

#ifndef REINISOLID_RK3_H_
#define REINISOLID_RK3_H_

class reinisolid_RK3 : public reinitopo, public increment
{
public:
	reinisolid_RK3(lexer* p);
	virtual ~reinisolid_RK3();
	virtual void start(lexer*,fdm*,ghostcell*,field&);

	int *sizeM;
	vec f,frk1,frk2,L,dt;

private:
	reinidisc *prdisc;

	void step(lexer*, fdm*);
    void time_preproc(lexer*);
	
	double starttime,endtime;

	int gcval,gcval_topo,gcval_initopo,reiniter,n;
	const double epsi;
};

#endif
