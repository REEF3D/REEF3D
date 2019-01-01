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
class fdm2D;
class ghostcell;
class slice;
class vec2D;
class cpt2D;

using namespace std;

#ifndef SOLVER2D_H_
#define SOLVER2D_H_

class solver2D
{

public:

	virtual void start(lexer*,fdm2D*,ghostcell*, slice&, vec2D&, vec2D&, int, int, double)=0;
	virtual void solve(lexer*,fdm2D*,ghostcell*, vec2D&, vec2D&, int, int, int&, int, double, cpt2D&)=0;
	virtual void setup(lexer*,fdm2D*,ghostcell*,int,cpt2D&)=0;	
};

#endif
