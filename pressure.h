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

class lexer;
class fdm;
class ghostcell;
class momentum;
class ioflow;
class poisson;
class solver;
class field;

using namespace std;

#ifndef PRESSURE_H_
#define PRESSURE_H_

class pressure
{
public:

	virtual void start(fdm*,lexer* p, poisson*, solver*, ghostcell*,momentum*,ioflow*,field&,field&,field&,double)=0;
	virtual void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double)=0;
	virtual void upgrad(lexer*,fdm*)=0;
	virtual void vpgrad(lexer*,fdm*)=0;
	virtual void wpgrad(lexer*,fdm*)=0;
	virtual void fillapu(lexer*,fdm*)=0;
	virtual void fillapv(lexer*,fdm*)=0;
	virtual void fillapw(lexer*,fdm*)=0;
	virtual void ptimesave(lexer*,fdm*,ghostcell*)=0;
};

#endif
