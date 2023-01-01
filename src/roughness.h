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

#include"increment.h"

class fdm;
class lexer;
class ghostcell;

using namespace std;

#ifndef ROUGHNESS_H_
#define ROUGHNESS_H_

class roughness : virtual public increment
{
public:
    roughness(lexer*);
	virtual ~roughness();

	virtual double ks_val(lexer*, fdm*,int,int,int,int,int);

private:
	double ks;

};

#endif


