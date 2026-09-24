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

#ifndef SFLOW_MOMENTUM_RK2_H_
#define SFLOW_MOMENTUM_RK2_H_

#include"sflow_momentum_func.h"
#include"slice4.h"

using namespace std;

class sflow_momentum_RK2 final : public sflow_momentum_func
{
public:
	sflow_momentum_RK2(lexer*, fdm2D*, ghostcell*, sflow_HLL*, sflow_signal_speed*, sflow_reconstruct*, sflow_diffusion*, sflow_pressure*, 
                        solver2D*, solver2D*, ioflow*, sflow_fsf*, sflow_forcing*, sixdof*);
	virtual ~sflow_momentum_RK2();
    
	void start(lexer*, fdm2D*, ghostcell*) override final;
    
private:
    slice4 WLRK1;
    slice4 UHRK1;
	slice4 VHRK1;
	slice4 WHRK1;
};

#endif
