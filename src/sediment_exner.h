/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include"vec2D.h"
#include"matrix2D.h"

class bedconc;
class topo_relax;
class solver2D;
class turbulence;
class ghostcell;

using namespace std;

#ifndef SEDIMENT_EXNER_H_
#define SEDIMENT_EXNER_H_

class sediment_exner : public topo, public increment
{
public:
	sediment_exner(lexer*, ghostcell*);
	virtual ~sediment_exner();
	virtual void start(lexer*, ghostcell*, sediment_fdm*);


private:
    void  topovel(lexer*,ghostcell*,sediment_fdm*,double&,double&,double&);
    void  timestep(lexer*,ghostcell*,sediment_fdm*);
    void  non_equillibrium_solve(lexer*,ghostcell*,sediment_fdm*);
    double  susp_qb(lexer*,ghostcell*,sediment_fdm*);
    
    bedconc *pcb;
    topo_relax *prelax;
    sediment_exnerdisc *pdx;
    solver2D *psolv;
    
    vec2D xvec,rhsvec;

	matrix2D M;
    
	int gcval_topo;
	double starttime;
    double maxdh,maxvz;
	double vx,vy,vz;
	double vzmax;
    double rhosed, rhowat, g, d50;
    double Ls;
    double tau_eff, shearvel_eff, shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
    
    slice4 q0,dqx0,dqy0;
};

#endif

