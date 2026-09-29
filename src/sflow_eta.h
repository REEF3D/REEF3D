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

#ifndef SFLOW_ETA_H_
#define SFLOW_ETA_H_

#include"sflow_fsf.h"
#include"increment.h"
#include"slice4.h" 

class patchBC_interface;

using namespace std;

class sflow_eta final : public sflow_fsf, public increment
{
public:    
	sflow_eta(lexer*, fdm2D*, ghostcell*,patchBC_interface*);
	virtual ~sflow_eta();
    
	void ini(lexer*, fdm2D*, ghostcell*, ioflow*) override final;
    void update(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&, slice&, double) override final;
	void depth_update(lexer*, fdm2D*, ghostcell*, slice&) override final;
	void wetdry(lexer*, fdm2D*, ghostcell*, slice&) override final;
    void wetdry_fluxes(lexer*, fdm2D*, ghostcell*, slice&) override final;
    void breaking(lexer*, fdm2D*, ghostcell*, slice&, slice&, double) override final;
    void breaking_persist(lexer*, fdm2D*, ghostcell*, slice&, slice&, double) override final;

private:
    void wetdrydeep(lexer*, fdm2D*, ghostcell*, slice&);
    
	int gcval_eta;
	double starttime;
    double wd_criterion;
    const double eps;
    
    patchBC_interface *pBC;
	slice4 K;
    int *temp;
};

#endif
