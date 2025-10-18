/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef NHFLOW_MOMENTUM_RK3_H_
#define NHFLOW_MOMENTUM_RK3_H_

#include"slice4.h"
#include"nhflow_momentum_func.h"
#include<vector>

class wind;
class vrans;

using namespace std;

class nhflow_momentum_RK3 : public nhflow_momentum_func
{
public:
	nhflow_momentum_RK3(lexer*, fdm_nhf*, ghostcell*, sixdof*, vrans*, nhflow_forcing*);
	virtual ~nhflow_momentum_RK3();
    
	void start(lexer*, fdm_nhf*, ghostcell*, ioflow*, nhflow_signal_speed*, nhflow_reconstruct*, nhflow_convection*, nhflow_diffusion*, nhflow_pressure*, solver*, solver*, nhflow*, nhflow_fsf*, nhflow_turbulence*,  vrans*) override;

    double *UHDIFF;
    double *VHDIFF;
    double *WHDIFF;
    
    double *UHRK1,*UHRK2;
    double *VHRK1,*VHRK2;
    double *WHRK1,*WHRK2;
    
    slice4 WLRK1,WLRK2;

private:
	
	int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
	double starttime;
    
    sixdof *p6dof;
    nhflow_forcing *pnhfdf;
    wind *pwind;
    vrans* pvrans;
};

#endif
