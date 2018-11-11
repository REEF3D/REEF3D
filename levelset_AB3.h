/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"freesurface.h"
#include"increment.h"
#include"field4.h"

class picard;
class fluid_update;
class heat;
class concentration;

using namespace std;

#ifndef LEVELSET_AB3_H_
#define LEVELSET_AB3_H_

class levelset_AB3 : public freesurface, increment
{
public:
	levelset_AB3(lexer*, fdm*, ghostcell*, heat*&, concentration*&);
	virtual ~levelset_AB3();
	virtual void start(fdm*,lexer*, discrete*, solver*, ghostcell*,ioflow*, reini*, particlecorr*,field&);
	virtual void ltimesave(lexer*,fdm*,field&);
    virtual void update(lexer*,fdm*,ghostcell*,field&);

	field4 lab1,lab2;

private:
    fluid_update *pupdate;
    picard *ppicard;

	int gcval_phi;
	double starttime;
};

#endif
