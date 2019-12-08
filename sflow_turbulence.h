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
Author: Hans Bihs
--------------------------------------------------------------------*/

class fdm2D;
class lexer;
class sflow_convection;
class sflow_diffusion;
class solver2D;
class ghostcell;
class ioflow;
class slice;

#ifndef SFLOW_TURBULENCE_H_
#define SFLOW_TURBULENCE_H_

using namespace std;

class sflow_turbulence
{

public:
	virtual void start(lexer*, fdm2D*, ghostcell*, sflow_convection*, sflow_diffusion*, solver2D*, ioflow*)=0;
	virtual void ktimesave(lexer*, fdm2D*, ghostcell*)=0;
	virtual void etimesave(lexer*, fdm2D*, ghostcell*)=0;
	
};

#endif
