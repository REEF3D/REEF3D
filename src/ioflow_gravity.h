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

#ifndef IOFLOW_GRAVITY_H_
#define IOFLOW_GRAVITY_H_

#include"ioflow.h"
#include"increment.h"

class fdm_fnpf;

using namespace std;


class ioflow_gravity final : public ioflow, public increment
{

public:

	ioflow_gravity(lexer*,ghostcell*,patchBC_interface*);
	virtual ~ioflow_gravity();
	void gcio_update(lexer*,fdm*,ghostcell*) override final;
    void gcio_update_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
	void inflow_walldist(lexer*,fdm*,ghostcell*,convection*,reini*,ioflow*) override final;
	void discharge(lexer*,fdm*,ghostcell*) override final;
	void inflow(lexer*,fdm*,ghostcell*,field&,field&,field&) override final;
	void rkinflow(lexer*,fdm*,ghostcell*,field&,field&,field&) override final;
	void fsfinflow(lexer*,fdm*,ghostcell*) override final;
	void fsfrkin(lexer*,fdm*,ghostcell*,field&) override final;
	void fsfrkout(lexer*,fdm*,ghostcell*,field&) override final;
    void iogcb_update(lexer*,fdm*,ghostcell*) override final;
    void isource(lexer*,fdm*,ghostcell*,vrans*) override final;
    void jsource(lexer*,fdm*,ghostcell*,vrans*) override final;
    void ksource(lexer*,fdm*,ghostcell*,vrans*) override final;
    void pressure_io(lexer*,fdm*,ghostcell*) override final;
    void turbulence_io(lexer*,fdm*,ghostcell*) override final;
    void veltimesave(lexer*,fdm*,ghostcell*,vrans*) override final;
    void flowfile(lexer*,fdm*,ghostcell*,turbulence*) override final;
    
    void wavegen_precalc(lexer*,ghostcell*) override final;
    void wavegen_precalc_ini(lexer*,ghostcell*) override final;
    void u_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void v_relax(lexer*,fdm*, ghostcell*,field&) override final;
    void w_relax(lexer*,fdm*, ghostcell*,field&) override final;
    void p_relax(lexer*,fdm*,ghostcell*,field&) override final;
	void phi_relax(lexer*,ghostcell*,field&) override final;
    void vof_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void turb_relax(lexer*,fdm*,ghostcell*,field&) override final;
    void U_relax(lexer*,ghostcell*,double*,double*) override final;
    void V_relax(lexer*,ghostcell*,double*,double*) override final;
    void W_relax(lexer*,ghostcell*,double*,double*) override final;
    void P_relax(lexer*,ghostcell*,double*) override final;
    void WL_relax(lexer*,ghostcell*,slice&,slice&) override final;
    void fi_relax(lexer*,ghostcell*,field&,field&) override final;
    void fivec_relax(lexer*, ghostcell*, double*) override final;
    void fifsf_relax(lexer*, ghostcell*, slice&) override final;
    void visc_relax(lexer*, ghostcell*, slice&) override final;
    void eta_relax(lexer*,ghostcell*,slice&) override final;
    void um_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
    void vm_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
	void wm_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
    void ws_relax(lexer*,ghostcell*,slice&,slice&,slice&) override final;
	void pm_relax(lexer*,ghostcell*,slice&) override final;
    
    void wavegen_2D_precalc(lexer*,fdm2D*,ghostcell*) override final;
    void wavegen_2D_precalc_ini(lexer*,ghostcell*) override final;
    
    void discharge2D(lexer*,fdm2D*,ghostcell*) override final;
    void waterlevel2D(lexer*,fdm2D*,ghostcell*,slice&) override final;
    void Qin2D(lexer*,fdm2D*,ghostcell*) override final;
	void Qout2D(lexer*,fdm2D*,ghostcell*) override final;
    void inflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&) override final;
	void rkinflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&) override final;
	void isource2D(lexer*,fdm2D*,ghostcell*) override final;
    void jsource2D(lexer*,fdm2D*,ghostcell*) override final;
	void full_initialize2D(lexer*,fdm2D*,ghostcell*) override final;
    
    double wave_fsf(lexer*,ghostcell*,double) override final;
    double wave_xvel(lexer*,ghostcell*,double,double,double) override final;
    double wave_yvel(lexer*,ghostcell*,double,double,double) override final;
    double wave_zvel(lexer*,ghostcell*,double,double,double) override final;
    
	int iozonecheck(lexer*,fdm*) override final;
    
    void ini(lexer*,fdm*,ghostcell*) override final;
    
    void waterlevel_update(lexer*,fdm*,ghostcell*) override final {};
    
    
    // fnpf
    void wavegen_precalc_fnpf(lexer*,fdm_fnpf*,ghostcell*) override final {};
    void ini_fnpf(lexer*,fdm_fnpf*,ghostcell*) override final;
    void inflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,double*,double*,slice&,slice&) override final;
    void rkinflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&) override final;
    void ini2D(lexer*,fdm2D*,ghostcell*) override final;
    void ini_ptf(lexer*,fdm*,ghostcell*) override final;
    
    // nhflow
    void wavegen_precalc_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void wavegen_precalc_ini_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void ini_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void discharge_nhflow(lexer*,fdm_nhf*,ghostcell*) override final;
    void inflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override final;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override final;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*) override final {};
    void isource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans_nhflow*,slice&) override final;
    void jsource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans_nhflow*,slice&) override final;
    void ksource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans_nhflow*,slice&) override final;
    void fsfinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,slice&) override final;
    void turb_relax_nhflow(lexer*,fdm_nhf*,ghostcell*,double*) override final {};
    
    void vrans_sed_update(lexer*,fdm*,ghostcell*,vrans*) override final;
	
	
private:
	double omega_x,omega_y;
	double theta_x,theta_y;
	double dist_x,dist_y,dist_z;
    
    patchBC_interface *pBC;
};

#endif

