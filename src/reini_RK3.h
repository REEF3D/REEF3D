/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef REINI_RK3_H_
#define REINI_RK3_H_

#include"reini.h"
#include"field4.h"
#include"increment.h"

class reinidisc;
class picard;

using namespace std;

class reini_RK3 : public reini, public increment
{
public:
	reini_RK3(lexer* p,int);
	virtual ~reini_RK3();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);

	int *sizeM;
	field4 frk1,frk2,dt;

private:
    picard *ppicard;
	reinidisc *prdisc;

	void step(lexer*, fdm*);
    void time_preproc(lexer*);

	
	double starttime,endtime;

	int gcval_phi,gcval_ro,gcval_iniphi,reiniter,n, gcval;
	const double epsi;
};

#endif
