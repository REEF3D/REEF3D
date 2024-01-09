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

#include"freesurface.h"
#include"gradient.h"

class picard;
class heat;
class concentration;
class fluid_update;

using namespace std;

#ifndef LEVELSET_RK2_H_
#define LEVELSET_RK2_H_

class levelset_RK2 : public freesurface, gradient
{
public:
	levelset_RK2(lexer*, fdm*, ghostcell*, heat*&, concentration*&);
	virtual ~levelset_RK2();
	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
    virtual void update(lexer*,fdm*,ghostcell*,field&);

private:
    fluid_update *pupdate;
    picard *ppicard;

	int gcval_phi;
	double starttime;
};
#endif
