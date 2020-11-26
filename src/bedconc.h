/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
#include"bedshear.h"

class lexer;
class fdm;
class sl;
class ghostcell;

#ifndef BEDCONC_H_
#define BEDCONC_H_

using namespace std;

class bedconc : public bedshear
{
public:
	bedconc(lexer*,turbulence*);
	virtual ~bedconc();
	double cbed(lexer*,fdm*,ghostcell*,field&);


private:
	int ii,jj,kk;
	int count,q;
	double ws,d50,ks,shields,kappa;
	double Rstar, g, visc;
	double rhosed,rhowat;
	double tau, taucrit;
	double val;
	double Ti,Ds;
    double tau_eff, shearvel_eff, shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
};
#endif



