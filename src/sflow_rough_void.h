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

#ifndef SFLOW_ROUGH_VOID_H_
#define SFLOW_ROUGH_VOID_H_

#include"sflow_roughness.h"
#include"increment.h"
#include"slice4.h"

using namespace std;

class sflow_rough_void : public sflow_roughness
{

public:
    sflow_rough_void(lexer*);
	virtual ~sflow_rough_void();
    
	virtual void u_source(lexer*, fdm2D*, slice&);
    virtual void v_source(lexer*, fdm2D*, slice&);

};

#endif
