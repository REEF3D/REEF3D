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

#include"nhflow_reconstruct.h"
#include"weno_nug_func.h"
#include"slice4.h"

class lexer;
class ghostcell;
class fdm_nhf;
class slice;
class patchBC_interface;

#ifndef NHFLOW_RECONSTRUCT_WENO_H_
#define NHFLOW_RECONSTRUCT_WENO_H_

using namespace std;

class nhflow_reconstruct_weno : public nhflow_reconstruct, public weno_nug_func
{
public:
	nhflow_reconstruct_weno(lexer*,patchBC_interface*);
	virtual ~nhflow_reconstruct_weno();

    virtual void reconstruct_2D_x(lexer*,ghostcell*,fdm_nhf*,slice&,slice&,slice&);
    virtual void reconstruct_2D_y(lexer*,ghostcell*,fdm_nhf*,slice&,slice&,slice&);
    virtual void reconstruct_2D_WL(lexer*,ghostcell*,fdm_nhf*);
    
    virtual void reconstruct_3D_x(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*);
    virtual void reconstruct_3D_y(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*);
    virtual void reconstruct_3D_z(lexer*,ghostcell*,fdm_nhf*,double*,double*,double*);
    
    slice4 dfdx;

private:
    void iqmin(lexer*, double*);
	void jqmin(lexer*, double*);
	void kqmin(lexer*, double*);
	void iqmax(lexer*, double*);
	void jqmax(lexer*, double*);
	void kqmax(lexer*, double*);
    
    void iqmin_sl(lexer*, slice&);
    void iqmax_sl(lexer*, slice&);
    void jqmin_sl(lexer*, slice&);
    void jqmax_sl(lexer*, slice&);
    
    double limiter(double, double);

    double ivel1,ivel2,jvel1,jvel2;
    double val,denom;
    double dfdx_min, dfdx_plus, dfdy_min, dfdy_plus, dfdz_min, dfdz_plus;
    int qq;
    double r, phi,minphi;
    
    double *DFDX;
    
    patchBC_interface *pBC;
};

#endif
