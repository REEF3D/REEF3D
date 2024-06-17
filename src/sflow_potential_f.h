/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef SFLOW_POTENTIAL_F_H_
#define SFLOW_POTENTIAL_F_H_

#include"sflow_potential.h"
#include"increment.h"
#include"sliceint4.h"

class slice;

using namespace std;

class sflow_potential_f : public sflow_potential, public increment
{

public:
	sflow_potential_f(lexer*);
	virtual ~sflow_potential_f();

	virtual void start(lexer*,fdm2D*, solver2D*, ghostcell*);


private:
	void ucalc(lexer*,fdm2D*,slice&);
	void vcalc(lexer*,fdm2D*,slice&);
    
    void laplace(lexer*,fdm2D*,slice&);
    void ini_bc(lexer*,fdm2D*,ghostcell*);
    
    
	double starttime,endtime;
	int count;
	int gcval_pot;
    
    sliceint4 bc;
};

#endif

