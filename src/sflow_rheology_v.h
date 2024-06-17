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

#ifndef SFLOW_RHEOLOGY_V_H_
#define SFLOW_RHEOLOGY_V_H_

#include"sflow_rheology.h"
#include"increment.h"
#include"slice4.h"

using namespace std;

class sflow_rheology_v : public sflow_rheology, public increment
{

public:
    sflow_rheology_v(lexer*);
	virtual ~sflow_rheology_v();
    
	virtual void u_source(lexer*, fdm2D*, slice&, slice&);
    virtual void v_source(lexer*, fdm2D*, slice&, slice&);

private:
    

};

#endif
