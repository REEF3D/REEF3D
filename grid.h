/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;

#ifndef GRID_H_
#define GRID_H_

using namespace std;

class grid :  public increment
{
public:

	grid (lexer *);
	virtual ~grid();
	
	void makegrid(lexer*,fdm*,ghostcell*);
	void update_topo_grid(lexer*,fdm*,ghostcell*);
	void update_sixdof_grid(lexer*,fdm*,ghostcell*);

   
	
};

#endif






