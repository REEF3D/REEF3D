/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"increment.h"

class lexer;
class fdm;
class field;

#ifndef SUSPRHS_H_
#define SUSPRHS_H_

using namespace std;

class susprhs : public increment
{
public:
	susprhs(lexer*);
	virtual ~susprhs();
	void suspsource(lexer*,fdm*,field&);
	void sedfsf(lexer*,fdm*,field&);
	void clearrhs(lexer*,fdm*);

private:
	int ii,jj,kk;
	int count,q;
	double ws,d50,ks,gi;
	double rhosed,rhowat;
};
#endif


