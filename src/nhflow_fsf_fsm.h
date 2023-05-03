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

#include"nhflow_fsf.h"
#include"increment.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class picard;
class fluid_update;
class heat;
class concentration;
class sflow_eta_disc;
class sflow_hxy_disc;
class patchBC_interface;

using namespace std;

#ifndef NHFLOW_FSF_FSM_H_
#define NHFLOW_FSF_FSM_H_

class nhflow_fsf_fsm : public nhflow_fsf, public increment
{
public:
    nhflow_fsf_fsm(lexer*, fdm_nhf*, ghostcell*,ioflow*,patchBC_interface*);
	virtual ~nhflow_fsf_fsm();
    
    virtual void start(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    
    virtual void step1(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    virtual void step2(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    virtual void step3(lexer*, fdm_nhf*, ghostcell*, ioflow*, double*, double*, double*, slice&, slice&, double);
    
    void update(lexer*,fdm_nhf*,ghostcell*,slice&);
    
private: 
    fluid_update *pupdate;
    sflow_eta_disc *peta;
	sflow_hxy_disc *phxy;
    patchBC_interface *pBC;
    
    slice1 P;
    slice2 Q;    
    
    int gcval_phi;
	double starttime;
    double phival,H;
	double d;
    const double epsi;

};

#endif
