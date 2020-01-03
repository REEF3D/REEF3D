/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
#include"ddweno_nug.h"
#include"vec.h"

class picard;
class heat;
class concentration;
class fluid_update;
class reinidisc;
class flux;

using namespace std;

#ifndef LEVELSET_RK3_V_H_
#define LEVELSET_RK3_V_H_

class levelset_RK3_V : public freesurface, ddweno_nug
{
public:
	levelset_RK3_V(lexer*, fdm*, ghostcell*, heat*&, concentration*&);
	virtual ~levelset_RK3_V();
	virtual void start(fdm*,lexer*, convection*, solver*, ghostcell*,ioflow*, reini*, particlecorr*,field&);
	virtual void ltimesave(lexer*,fdm*,field&);
    virtual void update(lexer*,fdm*,ghostcell*,field&);

private:
    void disc(lexer*,fdm*,vec&);
    
    fluid_update *pupdate;
    picard *ppicard;
    reinidisc *prdisc;
    flux *pflux;
    
    vec ark1,ark2,f,L;

	int gcval_phi;
	double starttime;
};
#endif
