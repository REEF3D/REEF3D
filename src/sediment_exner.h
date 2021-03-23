/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"topo.h"
#include"slice4.h"

class bedconc;
class topo_relax;
class turbulence;
class ghostcell;
class sediment_exnerdisc;

using namespace std;

#ifndef SEDIMENT_EXNER_H_
#define SEDIMENT_EXNER_H_

class sediment_exner : public topo, public increment
{
public:
	sediment_exner(lexer*, fdm*, ghostcell*,turbulence*);
	virtual ~sediment_exner();
	virtual void start(fdm*,lexer*, convection*, ghostcell*,reinitopo*);


private:
    void  topovel(lexer*,fdm*,ghostcell*,double&,double&,double&);
    void  timestep(lexer*,fdm*,ghostcell*);
    void  non_equillibrium_ini(lexer*,fdm*,ghostcell*);
    void  non_equillibrium_solve(lexer*,fdm*,ghostcell*);
    
    bedconc *pcb;
    topo_relax *prelax;
    sediment_exnerdisc *pdx;
    
	int gcval_topo;
	double starttime;
    double maxdh;
	double vx,vy,vz;
	double vzmax;
    double ws;
    double rhosed, rhowat, g, d50;
    
    slice4 dh;
    slice4 q0,dqx0,dqy0;
};

#endif

