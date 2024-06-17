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

#ifndef SFLOW_PRESSURE_H_
#define SFLOW_PRESSURE_H_

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;
class ioflow;
class patchBC_interface;

using namespace std;

class sflow_pressure
{
public:

	virtual void start(lexer*, fdm2D*, ghostcell*, solver2D*, ioflow*, slice&, slice&, slice&, slice&, slice&, slice&, double)=0;
	virtual void upgrad(lexer*, fdm2D*, slice&, slice&)=0;
	virtual void vpgrad(lexer*, fdm2D*, slice&, slice&)=0;
    virtual void wpgrad(lexer*, fdm2D*, slice&, slice&)=0;
    
    virtual void ucorr(lexer*,fdm2D*,slice&,slice&,double)=0;
	virtual void vcorr(lexer*,fdm2D*,slice&,slice&,double)=0;
	virtual void wcorr(lexer*,fdm2D*,double,slice&,slice&,slice&)=0;
    virtual void wcalc(lexer*,fdm2D*,double,slice&,slice&,slice&)=0;
};

#endif
