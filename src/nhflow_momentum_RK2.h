/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef NHFLOW_MOMENTUM_RK2_H_
#define NHFLOW_MOMENTUM_RK2_H_

#include"slice4.h"
#include"nhflow_momentum_func.h"
#include"nhflow_breaking.h"
#include<vector>

class wind;
class vrans;
class sediment;

using namespace std;

class nhflow_momentum_RK2 final : public nhflow_momentum_func, public nhflow_breaking
{
public:
	nhflow_momentum_RK2(lexer*, fdm_nhf*, ghostcell*, sixdof*, vrans_nhflow*, nhflow_forcing*, sediment*);
	virtual ~nhflow_momentum_RK2();
    
	void start(lexer*, fdm_nhf*, ghostcell*, ioflow*, nhflow_signal_speed*, nhflow_reconstruct*, nhflow_convection*, 
                nhflow_diffusion*, nhflow_pressure*, solver*, solver*, nhflow*, nhflow_fsf*, nhflow_turbulence*, vrans_nhflow*) override final;
    
    int stages() const override final { return 2; }
    void step_begin(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&) override final;
    void phase_F(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_M(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_P1(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    void phase_P2(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    double stage_alpha(int s) const override final { return s==0?1.0:0.5; }
    void phase_E(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) override final;
    slice& stage_WL(fdm_nhf*,int) override final;
    double* stage_UH(fdm_nhf*,int,int) override final;

    double *UHDIFF,*VHDIFF,*WHDIFF;
    double *UHRK1,*VHRK1,*WHRK1;
    
    slice4 WLRK1;

private:
	int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
	double starttime;
    
    nhflow_convection *pweno;
    sixdof *p6dof;
    nhflow_forcing *pnhfdf;
    wind *pwind;
    vrans_nhflow* pvrans;
    sediment *psed;
};

#endif
