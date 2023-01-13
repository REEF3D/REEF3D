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

#include"reini.h"
#include"ddweno.h"
#include"vec.h"
#include"increment.h"

class reinidisc;
class picard;

using namespace std;

#ifndef REINI_RK3_H_
#define REINI_RK3_H_

class reini_RK3 : public reini, public increment
{
public:
	reini_RK3(lexer* p,int);
	virtual ~reini_RK3();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);
    virtual void startV(fdm*,lexer*,vec&,ghostcell*,ioflow*);

	int *sizeM;
	vec f,frk1,frk2,L,dt;

private:
    picard *ppicard;
	reinidisc *prdisc;

    void fsfrkioV(lexer*, fdm*, ghostcell*,vec&);
	void step(lexer*, fdm*);
    void time_preproc(lexer*);
	void inisolid(lexer*, fdm*);

	
	double starttime,endtime;

	int gcval_phi,gcval_ro,gcval_iniphi,reiniter,n;
	const double epsi;
};

#endif
