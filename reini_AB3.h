/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"reini.h"
#include"vec.h"
#include"increment.h"

class reinidisc;
class picard;

using namespace std;

#ifndef REINI_AB3_H_
#define REINI_AB3_H_

class reini_AB3 : public reini, increment
{
public:
	reini_AB3(lexer* p, fdm *a);
	virtual ~reini_AB3();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);
    virtual void startV(fdm*,lexer*,vec&,ghostcell*,ioflow*);

	int *sizeM;
    vec f,dab1,dab2,L,dt;


private:
    picard *ppicard;
	reinidisc *prdisc;

	void fsfrkioV(lexer*, fdm*, ghostcell*,vec&);
	void step(lexer*, fdm*);
    void time_preproc(lexer*);
	void inisolid(lexer*, fdm*);
	double minsign,maxdiff;

	int gcval_phi,gcval_iniphi,gcval_ro,reiniter;	
	double starttime;
};

#endif
