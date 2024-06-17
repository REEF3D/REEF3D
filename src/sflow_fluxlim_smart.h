/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef SFLOW_FLUXLIM_SMART_H_
#define SFLOW_FLUXLIM_SMART_H_

#include"sflow_fluxlim.h"
#include"increment.h"

class slice;
class fdm;
class lexer;

using namespace std;

class sflow_fluxlim_smart : public sflow_fluxlim, public increment 
{
public:
    sflow_fluxlim_smart (lexer *);
	virtual ~sflow_fluxlim_smart();

	virtual double iphi(slice&,int,int,int,int);
	virtual double jphi(slice&,int,int,int,int);
    
private:
    double r, phi,denom;
	double dx,dy;
    double minphi;
	double L;

};

#endif

