/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef VRANS_NHFLOW_BASE_H_
#define VRANS_NHFLOW_BASE_H_

#include"vrans.h"
#include"increment.h"

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

class vrans_nhflow_base
{
public:
	virtual void initialize(lexer*, fdm_nhf*, ghostcell*)=0;	
	virtual void start(lexer*, fdm_nhf*, ghostcell*, int)=0;

	virtual void u_source(lexer*, fdm_nhf*)=0;
	virtual void v_source(lexer*, fdm_nhf*)=0;
	virtual void w_source(lexer*, fdm_nhf*)=0;
    
    virtual void ke_source(lexer*, fdm_nhf*, field&)=0;
    virtual void kw_source(lexer*, fdm_nhf*, field&)=0;
    virtual void eps_source(lexer*, fdm_nhf*, field&, field&)=0;
    virtual void omega_source(lexer*, fdm_nhf*, field&, field&)=0;
    
    virtual void eddyv_func(lexer*, fdm_nhf*)=0;
    

};

#endif
