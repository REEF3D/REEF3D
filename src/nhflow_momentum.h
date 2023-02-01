/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

class lexer;
class fdm_nhf;
class ghostcell;
class nhflow_convection;
class diffusion;
class nhflow_pressure;
class turbulence;
class solver;
class ioflow;
class nhflow;
class nhflow_fsf;
class vrans;

using namespace std;

#ifndef NHFLOW_MOMENTUM_H_
#define NHFLOW_MOMENTUM_H_

class nhflow_momentum
{
public:

	virtual void start(lexer*, fdm_nhf*, ghostcell*, ioflow*, nhflow_convection*, diffusion*, nhflow_pressure*, solver*, nhflow*, nhflow_fsf*, vrans*)=0;

};

#endif
