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

#include"benchmark.h"
#include"increment.h"
#include"gradient.h"
#include"field4.h"

class fdm;
class lexer;
class ghostcell;

#ifndef BENCHMARK_TAYLORGREEN_H_
#define BENCHMARK_TAYLORGREEN_H_

using namespace std;

class benchmark_TaylorGreen : public benchmark, public increment, public gradient
{

public:
    benchmark_TaylorGreen(lexer*,fdm*);
	virtual ~benchmark_TaylorGreen();

	virtual void start(lexer*, fdm*, ghostcell*, convection*);
private:
	
	field4 vx,vy,vz;
};

#endif




