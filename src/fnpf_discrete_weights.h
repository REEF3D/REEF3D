/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

/*
    Reference:
    Bengt Fornberg, Generation of Finite Difference Formulas on Arbitrarily Spaced Grids, 
    Mathematics of Compuation, 51, 184, 1988, pp. 699-706
*/

#include"fnpf_convection.h"
#include"increment.h"

#ifndef FNPF_DISCRETE_WEIGHTS_H_
#define FNPF_DISCRETE_WEIGHTS_H_

using namespace std;

class fnpf_discrete_weights : public increment
{
public:
	fnpf_discrete_weights(lexer*);
	virtual ~fnpf_discrete_weights();

    void ck_weights(lexer*, double**, double *, int, int, int, int);

private:
    

};

#endif









