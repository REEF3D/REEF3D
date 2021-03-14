/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"bedconc.h"

class lexer;
class fdm;
class field;
class turbulence;

#ifndef IBCSUSP_H_
#define IBCSUSP_H_

using namespace std;

class ibcsusp : public bedconc
{
public:
	ibcsusp(lexer*,turbulence*);
	virtual ~ibcsusp();
	void ibcsusp_start(lexer*,fdm*,ghostcell*,field&);


private:
	int ii,jj,kk;
	int count,q,n;
	const double epsi;
	double g,d50;
	double shields,Rstar,visc;
	double rhosed,rhowat;
    double kappa,u_plus,u_abs;
    double tau,taucrit,ks;

};
#endif
