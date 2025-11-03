/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef SFLOW_MOMENTUM_FUNC_H_
#define SFLOW_MOMENTUM_FUNC_H_

#include"sflow_momentum.h"
#include"sflow_bcmom.h"

class sflow_fsf;
class sflow_signal_speed;
class sflow_reconstruct;
class sflow_fsf_reconstruct;

using namespace std;

class sflow_momentum_func : public sflow_momentum, public sflow_bcmom
{
public:
	sflow_momentum_func(lexer*, fdm2D*, ghostcell*);
	virtual ~sflow_momentum_func();
    

    virtual void inidisc(lexer*, fdm2D*, ghostcell*, sflow_fsf*);
    void reconstruct(lexer*, fdm2D*, ghostcell*, sflow_fsf*, sflow_signal_speed*, sflow_reconstruct*,slice&,slice&,slice&,slice&,slice&,slice&,slice&);
    void velcalc(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&);
    
	void irhs(lexer*,fdm2D*,ghostcell*);
	void jrhs(lexer*,fdm2D*,ghostcell*);
	void krhs(lexer*,fdm2D*,ghostcell*);
    
    void clearrhs(lexer*,fdm2D*,ghostcell*);
	
    

	int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
    
    double starttime;

};

#endif
