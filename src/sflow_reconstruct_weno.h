/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef SFLOW_RECONSTRUCT_WENO_H_
#define SFLOW_RECONSTRUCT_WENO_H_

#include"sflow_reconstruct.h"
#include"weno_nug_func.h"
#include"slice4.h"

class lexer;
class ghostcell;
class fdm2D;
class slice;
class patchBC_interface;

using namespace std;

class sflow_reconstruct_weno : public sflow_reconstruct, public weno_nug_func
{
public:
	sflow_reconstruct_weno(lexer*,patchBC_interface*);
	virtual ~sflow_reconstruct_weno();

    virtual void reconstruct_x(lexer*,ghostcell*,fdm2D*,slice&,slice&,slice&);
    virtual void reconstruct_y(lexer*,ghostcell*,fdm2D*,slice&,slice&,slice&);
    virtual void reconstruct_WL(lexer*,ghostcell*,fdm2D*);
    
    slice4 dfdx,dfdy;

private:
    void iqmin(lexer*, slice&);
    void iqmax(lexer*, slice&);
    void jqmin(lexer*, slice&);
    void jqmax(lexer*, slice&);
    
    double limiter(double, double);

    double ivel1,ivel2,jvel1,jvel2;
    int qq;
    
    patchBC_interface *pBC;
};

#endif
