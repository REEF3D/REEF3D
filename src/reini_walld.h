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

#ifndef REINI_WALLD_H_
#define REINI_WALLD_H_

#include"reini.h"
#include"field4.h"
#include"gradient.h"

class reinidisc;
class picard;

using namespace std;

class reini_walld : public reini, gradient
{
public:
	reini_walld(lexer* p, fdm *a);
	virtual ~reini_walld();
	virtual void start(fdm*,lexer*,field&,ghostcell*,ioflow*);

    field4 dab;
	reinidisc *prdisc;

private:


	void step(fdm*, lexer*);
	double dt,minsign,maxdiff;

	int gcval_phi,gcval_iniphi,gcval_ro,reiniter;	
	int q;
	double starttime;
};

#endif
