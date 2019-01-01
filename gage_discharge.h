/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

#ifndef GAGE_DISCHARGE_H_
#define GAGE_DISCHARGE_H_

class gage_discharge : public boundarycheck
{
public:
    gage_discharge(lexer*,fdm*,ghostcell*);
	virtual ~gage_discharge();

	void start(lexer*, fdm*, ghostcell*);

private:
    void ini_location(lexer*, fdm*, ghostcell*);
    int conv(double);

    int *iloc,*flag;
    double *q;
	double area,Ai;
    int n;
    ofstream qout;

};

#endif
