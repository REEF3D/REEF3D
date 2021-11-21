/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

using namespace std;

#ifndef NHFLOW_FSF_F_H_
#define NHFLOW_FSF_F_H_

class nhflow_fsf_f : public nhflow_fsf, public increment
{
public:
    nhflow_fsf_f(lexer*, fdm*, ghostcell*,ioflow*);
	virtual ~nhflow_fsf_f();
    
    virtual void start(lexer*, fdm*, ghostcell*, ioflow*);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*);
    
	void ltimesave(lexer*,fdm*,slice&);
    void update(lexer*,fdm*,ghostcell*,slice&);
    
private: 
    fluid_update *pupdate;
    picard *ppicard;
    
    int gcval_phi;
	double starttime;
    double phival,H;
	double d;
    const double epsi;
	
	sflow_eta_disc *peta;
	sflow_hxy_disc *phxy;
};

#endif
