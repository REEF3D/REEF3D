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

#include"boundarycheck.h"
#include"bedshear.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;
class sediment;

using namespace std;

#ifndef BEDSHEAR_PROBE_H_
#define BEDSHEAR_PROBE_H_

class bedshear_probe : public boundarycheck
{
public:
    bedshear_probe(lexer*,fdm*,ghostcell*);
	virtual ~bedshear_probe();

	void bedshear_gauge(lexer*, fdm*, ghostcell*, sediment*);


private:
    void ini_location(lexer*, fdm*, ghostcell*);
    void write(lexer*, fdm*, ghostcell*);
    int conv(double);

    int *iloc,*jloc,*flag;
    double *bsg;
    int n;
    ofstream bsgout;

};

#endif
