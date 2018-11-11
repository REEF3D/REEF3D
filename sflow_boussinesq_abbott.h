/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sflow_boussinesq.h"
#include"increment.h"
#include"slice1.h"
#include"slice2.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;

using namespace std;

#ifndef SFLOW_BOUSSINESQ_ABBOTT_H_
#define SFLOW_BOUSSINESQ_ABBOTT_H_

class sflow_boussinesq_abbott : public sflow_boussinesq, public increment
{
public:

    sflow_boussinesq_abbott(lexer*, fdm2D*);
	~sflow_boussinesq_abbott();
    

	virtual void ini(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    
	virtual void psi1(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
	virtual void psi2(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
    

private:
    void psi1_calc(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
	void psi2_calc(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
    
    void psi1_update(lexer*, fdm2D*);
	void psi2_update(lexer*, fdm2D*);
    
    
    slice1 Pxx;
    slice1 Pxx_n;
    
    slice2 Qyy;
    slice2 Qyy_n;

    
    double ddx,ddy,d;
};

#endif
