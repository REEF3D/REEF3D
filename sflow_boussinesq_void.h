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
#include"slice1.h"
#include"slice2.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;

using namespace std;

#ifndef SFLOW_BOUSSINESQ_VOID_H_
#define SFLOW_BOUSSINESQ_VOID_H_

class sflow_boussinesq_void : public sflow_boussinesq
{
public:

    sflow_boussinesq_void(lexer*, fdm2D*);
	virtual ~sflow_boussinesq_void();

	virtual void ini(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    
	virtual void psi1(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
	virtual void psi2(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, double);
    
    

};

#endif
