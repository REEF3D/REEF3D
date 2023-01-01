/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"sflow_fsf.h"
#include"increment.h"
#include"slice4.h" 

class sflow_eta_weno;
class sflow_hxy_disc;
class patchBC_interface;

using namespace std;

#ifndef SFLOW_ETA_H_
#define SFLOW_ETA_H_

class sflow_eta : public sflow_fsf, public increment
{
public:    
	sflow_eta(lexer*, fdm2D*, ghostcell*,patchBC_interface*);
	virtual ~sflow_eta();
	
    virtual void start(lexer*, fdm2D*, ghostcell*, ioflow*,slice&,slice&,double);
	virtual void ini(lexer*, fdm2D*, ghostcell*, ioflow*);
	virtual void depth_update(lexer*, fdm2D*, ghostcell*,slice&,slice&,slice&,slice&);
    virtual void breaking(lexer*, fdm2D*, ghostcell*, slice&, slice&, double);
    virtual void breaking_persist(lexer*, fdm2D*, ghostcell*, slice&, slice&, double);
	virtual void wetdry(lexer*, fdm2D*, ghostcell*,slice&,slice&,slice&);
private:
    
	
	
	int gcval_eta;
    
    double hxp,hxm,hyp,hym;
	double starttime;
	
    double wd_criterion;
    
	sflow_eta_weno *pconvec;
	sflow_hxy_disc *phxy;
    patchBC_interface *pBC;
	
	slice4 Lab;
    
        

};

#endif
