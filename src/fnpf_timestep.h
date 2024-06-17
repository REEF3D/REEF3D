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

#ifndef FNPF_TIMESTEP_H_
#define FNPF_TIMESTEP_H_

#include"timestep.h"
#include"increment.h"

class fdm_fnpf;
class lexer;
class ghostcell;

using namespace std;

class fnpf_timestep : public increment
{
public:
	fnpf_timestep(lexer*);
	virtual ~fnpf_timestep();
	virtual void start(fdm_fnpf*, lexer*,ghostcell*);
	virtual void ini(fdm_fnpf*, lexer*,ghostcell*);


private:
	double sqd,wallu,wallv,wallw;
	double cu,cv,cw;
	const double epsi;
	double isor,jsor,ksor;
	double irsm,jrsm,krsm;
    const double maxtimestep, c0_orig;
    double dx;

};

#endif
