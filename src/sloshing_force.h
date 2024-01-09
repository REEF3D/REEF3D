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

#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

#ifndef SLOSHING_FORCE_H_
#define SLOSHING_FORCE_H_

class sloshing_force : public increment
{
public:
    sloshing_force(lexer*,fdm*,ghostcell*);
	virtual ~sloshing_force();

	void start(lexer*, fdm*, ghostcell*);
    void force(lexer*, fdm*, ghostcell*);


private:
    void write(lexer*, fdm*, ghostcell*);
    
    double Fx_l,Fx_r,Fz,M;

    ofstream result;

};

#endif
