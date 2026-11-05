/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef PRINT_WSFLINE_Y_H_
#define PRINT_WSFLINE_Y_H_

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class wsfline_core;
class fdm;
class ghostcell;
class field;
class ioflow;

using namespace std;

class print_wsfline_y : public boundarycheck
{
public:
    print_wsfline_y(lexer*,fdm*,ghostcell*);
	virtual ~print_wsfline_y();

	void wsfline(lexer*, fdm*, ghostcell*,ioflow*);


private:
    void ini_location(lexer*, fdm*, ghostcell*);
    int *iloc;
    int n,q;
    ofstream wsfout;

	
    wsfline_core *pcore;

};

#endif

