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

#include"sflow_timestep.h"
#include"increment.h"

using namespace std;

#ifndef SFLOW_ETIMESTEP_H_
#define SFLOW_ETIMESTEP_H_

class sflow_etimestep : public sflow_timestep, public increment
{
public:
	sflow_etimestep(lexer*,fdm2D*);
	virtual ~sflow_etimestep();
	
    virtual void start(lexer*,fdm2D*,ghostcell*);
	virtual void ini(lexer*,fdm2D*,ghostcell*);
	
private:
	double cu,cv,velmax,wd_criterion;

};

#endif

