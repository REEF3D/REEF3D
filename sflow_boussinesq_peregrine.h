/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#ifndef SFLOW_BOUSSINESQ_PEREGRINE_H_
#define SFLOW_BOUSSINESQ_PEREGRINE_H_

class sflow_boussinesq_peregrine : public sflow_boussinesq, public increment
{
public:

    sflow_boussinesq_peregrine(lexer*, fdm2D*);
	~sflow_boussinesq_peregrine();
    

	virtual void ini(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    
	virtual void psi1(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
	virtual void psi2(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
    

private:
    void psi1_calc(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
	void psi2_calc(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
    
    void psi1_update(lexer*, fdm2D*);
	void psi2_update(lexer*, fdm2D*);
    
    
    slice1 Pxx,Qxy;
    slice1 Pxx_n,Qxy_n;
    slice1 P1x,Q1x,Q1y;
    slice1 P1x_n,Q1x_n,Q1y_n;
    
    slice2 Qyy,Pxy;
    slice2 Qyy_n,Pxy_n;
    slice1 P2x,P2y,Q2x,Q2y;
    slice1 P2x_n,P2y_n,Q2x_n,Q2y_n;

    
    double ddx,ddy,d;
};

#endif
