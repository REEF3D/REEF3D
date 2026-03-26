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

#ifndef SFLOW_HYDROSTATIC_H_
#define SFLOW_HYDROSTATIC_H_

#include"sflow_pressure.h"
#include"increment.h"
#include"slice4.h"

using namespace std;

class sflow_hydrostatic final : public sflow_pressure, public increment
{
public:
    sflow_hydrostatic(lexer*, fdm2D*,patchBC_interface*);
	virtual ~sflow_hydrostatic();
    
	void start(lexer*, fdm2D*, ghostcell*, solver2D*, ioflow*, slice&, slice&, slice&, slice&, slice&, slice&, double) override final;
	void upgrad(lexer*, fdm2D*, slice&, slice&) override final;
	void vpgrad(lexer*, fdm2D*, slice&, slice&) override final;
    void wpgrad(lexer*, fdm2D*, slice&, slice&) override final;
    
    void ucorr(lexer*,fdm2D*,slice&,slice&,double) override final;
	void vcorr(lexer*,fdm2D*,slice&,slice&,double) override final;
	void wcorr(lexer*,fdm2D*,double,slice&,slice&,slice&) override final;
    void wcalc(lexer*,fdm2D*,double,slice&,slice&,slice&) override final;
    
private:
    patchBC_interface *pBC;

};

#endif
