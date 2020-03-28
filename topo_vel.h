/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"ccipol.h"
#include"norm_vec.h"

class bedconc;
class topo_relax;
class turbulence;
class ghostcell;
class fnpf_convection;

using namespace std;

#ifndef TOPO_VEL_H_
#define TOPO_VEL_H_

class topo_vel : public ccipol, public norm_vec
{
public:
	topo_vel(lexer*,turbulence*);
	virtual ~topo_vel();

    void  topovel(lexer*,fdm*,ghostcell*,double&,double&,double&);
    
    void  topovel_xy_cds(lexer*,fdm*,ghostcell*,double&,double&,double&);
	
	void rkio_update(lexer*,fdm*,ghostcell*,field&);

	double xc,yc,zc;
	double nx,ny,nz,norm;
	double a1x,a1y,a1z,a2x,a2y,a2z;
	double b1x,b1y,b1z,b2x,b2y,b2z;
	double ascale, bscale;
	double da,db;
	double qx1,qx2,qy1,qy2;
	double dh;
	double signx,signy;
	double ws;
	double rhosed, rhowat, g, d50;

	bedconc *pcb;
    topo_relax *prelax;
    fnpf_convection *pdx;

	const double dx,epsi;
};

#endif
