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
#include"sflow_sigma.h"

class sflow_fsf;
class sflow_signal_speed;
class sflow_reconstruct;
class sflow_fsf_reconstruct;

using namespace std;

class sflow_momentum_func : public sflow_momentum, public sflow_bcmom, public sflow_sigma
{
public:
	sflow_momentum_func(lexer*, fdm_nhf*, ghostcell*);
	virtual ~sflow_momentum_func();
    

    virtual void inidisc(lexer*, fdm_nhf*, ghostcell*, sflow_fsf*);
    void reconstruct(lexer*, fdm_nhf*, ghostcell*, sflow_fsf*, sflow_signal_speed*, sflow_reconstruct*,slice&,double*,double*,double*,double*,double*,double*);
    void velcalc(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,slice&);
    
	void irhs(lexer*,fdm_nhf*,ghostcell*);
	void jrhs(lexer*,fdm_nhf*,ghostcell*);
	void krhs(lexer*,fdm_nhf*,ghostcell*);
    
    void clearrhs(lexer*,fdm_nhf*,ghostcell*);
	
    

	int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
    
    double starttime;

};

#endif
