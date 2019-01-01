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

#include"bedconc.h"

class lexer;
class fdm;
class field;
class turbulence;

#ifndef BCSUSP_H_
#define BCSUSP_H_

using namespace std;

class bcsusp : public bedconc
{
public:
	bcsusp(lexer*,turbulence*);
	virtual ~bcsusp();
	void bcsusp_start(lexer*,fdm*,ghostcell*,field&);

private:
	int ii,jj,kk;
	int count,q,n;
	const double epsi;
	double d50,ks,gi;
	double rhosed,rhowat,Rstar;
	double shields,visc;
    double kappa,u_plus,g,u_abs;
    double tau,taucrit;

};
#endif

