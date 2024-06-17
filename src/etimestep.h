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

#ifndef ETIMESTEP_H_
#define ETIMESTEP_H_

#include"timestep.h"
#include"increment.h"

class turbulence;

using namespace std;

class etimestep : public timestep, public increment
{
public:
	etimestep(lexer*);
	virtual ~etimestep();
	virtual void start(fdm*,lexer*,ghostcell*,turbulence*);
	virtual void ini(fdm*,lexer*,ghostcell*);


private:
	double max(double,double,double);
	double max(double,double);
	double min(double,double,double);
	double min(double,double);

	double visccrit,sqd,wallu,wallv,wallw;
	double uplus;
	double cu,cv,cw,ck,ce;
	double velmax;
	const double epsi;
	double isor,jsor,ksor;
	double irsm,jrsm,krsm;
	const double c0_orig;
    double dx,visc;


};

#endif
