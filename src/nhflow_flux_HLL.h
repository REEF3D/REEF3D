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

#include"nhflow_flux_fsf.h"
#include"nhflow_flux_weno.h"
#include"slice1.h"
#include"slice2.h"


#ifndef NHFLOW_FLUX_HLL_H_
#define NHFLOW_FLUX_HLL_H_

using namespace std;

class nhflow_flux_HLL : public nhflow_flux_fsf, public nhflow_flux_weno
{
public:
	nhflow_flux_HLL(lexer*,patchBC_interface*);
	virtual ~nhflow_flux_HLL();

	virtual void face_flux_2D(lexer*,fdm_nhf*,slice&,slice&,slice&,slice&,slice&);
    virtual void face_flux_3D(lexer*,ghostcell*,fdm_nhf*,slice&,double*,double*,double*,double*);
    
    double *Fs,*Fn,*Fe,*Fw;
    slice1 ETAs,ETAn;
    slice2 ETAe,ETAw;
    

private:
    

    double ivel1,ivel2,jvel1,jvel2;
    int qq;
    
    patchBC_interface *pBC;
};

#endif
