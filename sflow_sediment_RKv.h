/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"sflow_sediment_RK.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;

using namespace std;

#ifndef SFLOW_SEDIMENT_RKv_H_
#define SFLOW_SEDIMENT_RKv_H_

class sflow_sediment_RKv : public sflow_sediment_RK
{
public:
    sflow_sediment_RKv(lexer*, fdm2D*);
	virtual ~sflow_sediment_RKv();

    virtual void step1(lexer*, fdm2D*, ghostcell*, slice&, slice&, double);
    virtual void step2(lexer*, fdm2D*, ghostcell*, slice&, slice&, double);
    virtual void step3(lexer*, fdm2D*, ghostcell*, slice&, slice&, double);
    virtual void step4(lexer*, fdm2D*, ghostcell*, slice&, slice&, double);
    
private:
    
};

#endif
