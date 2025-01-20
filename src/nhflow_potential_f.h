/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_potential.h"
#include"increment.h"

using namespace std;

#ifndef NHFLOW_POTENTIAL_F_H_
#define NHFLOW_POTENTIAL_F_H_

class nhflow_potential_f : public nhflow_potential, public increment
{

public:
	nhflow_potential_f(lexer*);
	virtual ~nhflow_potential_f();

	virtual void start(lexer*,fdm_nhf*, solver*, ghostcell*);


private:
    void rhs(lexer*,fdm_nhf*);
	void ucalc(lexer*,fdm_nhf*);
	void vcalc(lexer*,fdm_nhf*);
	void wcalc(lexer*,fdm_nhf*);
    
    void laplace(lexer*,fdm_nhf*,ghostcell*);
    void ini_bc(lexer*,fdm_nhf*,ghostcell*);
    
    
	double starttime,endtime;
	int count;
	int gcval_pot;
    double sigxyz2;
    
    double *PSI;
    
    int *BC;

};

#endif

