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

#include"freesurface.h"
#include"gradient.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class picard;
class heat;
class convection;
class fluid_update;

using namespace std;

#ifndef VOF_AB_H_
#define VOF_AB_H_

class VOF_AB : public freesurface, gradient
{
public:
	VOF_AB(lexer*, fdm*, ghostcell*,heat*);
	virtual ~VOF_AB();
	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particle_corr*,field&);
	virtual void update(lexer*,fdm*,ghostcell*,field&);

	void compression(lexer*,fdm*,ghostcell*,convection*,field&,double);
	
private:
    fluid_update *pupdate;
	
	field1 uc;
	field2 vc;
	field3 wc;
    field4 F;
	field4 lab;

	int gcval_frac;
	double starttime;
	
	convection *ppconvec;
};
#endif

