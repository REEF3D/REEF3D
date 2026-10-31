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

#ifndef NHFLOW_SUSPENDED_IM1_H_
#define NHFLOW_SUSPENDED_IM1_H_

#include"nhflow_suspended.h"
#include"increment.h"
#include"slice4.h"

using namespace std;

class nhflow_suspended_IM1 final : public nhflow_suspended, public increment
{
public:
	nhflow_suspended_IM1(lexer *);
	virtual ~nhflow_suspended_IM1();
	void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, sediment_fdm*) override final;
	void ctimesave(lexer*, fdm_nhf*) override final;
    double* ckpt_wln() override final {return WLN;}
    int* ckpt_wlini() override final {return &wl_ini;}
    
    void suspsource(lexer*,fdm_nhf*,double*,sediment_fdm*);
    void bcsusp_start(lexer*,fdm_nhf*,ghostcell*,sediment_fdm*,double*);
	void clearrhs(lexer*,fdm_nhf*);
    void fillconc(lexer*,fdm_nhf*,ghostcell*,sediment_fdm*);
    void inflow_conc(lexer*,fdm_nhf*,sediment_fdm*);
    int cin_on = 0;
    
    void bed_rate(lexer*, ghostcell*, sediment_fdm*) override final;
    void bed_done(lexer*, sediment_fdm*) override final;

	int gcval_susp;

private:
    void timesource(lexer*, fdm_nhf*, double*);
    void drysave(lexer*, fdm_nhf*, sediment_fdm*);
    double starttime;
    void fill_wvel(lexer*,fdm_nhf*,ghostcell*,sediment_fdm*); 
    double *WVEL;
    double *WLN;    // water depth belonging to the stored concentration (end of the last solve)
    int wl_ini;
    
    // exchange with the bed since the last bed update: ED = sum ws (cbe - c_0) dt [m], its mean rate RT,
    // acc[0] the flow time, acc[1] the number of flow steps (checkpoints: raw entry of the grid)
    void exchange_sum(lexer*, fdm_nhf*, sediment_fdm*);
    slice4 ED,RT;
    double acc[2];
    unsigned long long ckid;
    double Rouse_formula(lexer*,fdm_nhf*,sediment_fdm*,double); 
    
    int count,q;
};

#endif
