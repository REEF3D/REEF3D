/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"sflow_hxy_disc.h"
#include"increment.h"

class sflow_flux;
class sflow_fluxlim;

#ifndef SFLOW_HXY_HIRES_H_
#define SFLOW_HXY_HIRES_H_

using namespace std;

class sflow_hxy_hires : public sflow_hxy_disc, public increment
{

public:

	sflow_hxy_hires(lexer*,patchBC_interface*,int);
	virtual ~sflow_hxy_hires();

	virtual void start(lexer*,slice&,slice&,slice&,int*,slice&,slice&,slice&);
	
private:
    double fx(lexer*, slice&, int, double);
	double fy(lexer*, slice&, int, double);
	
	double ul,ur,vl,vr,wl,wr;

	double dx,dy,dz;
	double L;
	
	sflow_fluxlim *plim;
    sflow_flux *pflux;
    
    double ivel1,ivel2,jvel1,jvel2;
    
    patchBC_interface *pBC;
};

#endif

