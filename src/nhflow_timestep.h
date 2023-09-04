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

#include"timestep.h"
#include"increment.h"

class fdm_nhf;
class lexer;
class ghostcell;

using namespace std;

#ifndef NHFLOW_TIMESTEP_H_
#define NHFLOW_TIMESTEP_H_

class nhflow_timestep : public increment
{
public:
	nhflow_timestep(lexer*);
	virtual ~nhflow_timestep();
	virtual void start(lexer*,fdm_nhf*,ghostcell*);
	virtual void ini(lexer*,fdm_nhf*,ghostcell*);


private:
	double wallu,wallv,wallw;
	double cu,cv,cw,co;
	const double epsi;
	double isor,jsor,ksor;
	double irsm,jrsm,krsm;
    const double maxtimestep, c0_orig;
    double dx;

};

#endif
