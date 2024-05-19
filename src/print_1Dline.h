/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

#ifndef PRINT_1DLINE_H_
#define PRINT_1DLINE_H_

class print_1Dline : public boundarycheck
{
public:
    print_1Dline(lexer*,fdm*,ghostcell*);
	virtual ~print_1Dline();

	void height_gauge(lexer*, fdm*, ghostcell*);

private:
    void ini_location(lexer*, fdm*, ghostcell*);
    void write(lexer*, fdm*, ghostcell*);
    int conv(double);
    int *iloc,*jloc,*flag;
    double *wsf;
    int n;
    ofstream wsfout;
};

#endif


