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

#include"ioflow.h"
#include"increment.h"

class fdm_fnpf;

using namespace std;

#ifndef IOFLOW_GRAVITY_H_
#define IOFLOW_GRAVITY_H_


class ioflow_gravity : public ioflow, public increment
{

public:

	ioflow_gravity(lexer*,ghostcell*,patchBC_interface*);
	virtual ~ioflow_gravity();
	virtual void gcio_update(lexer*,fdm*,ghostcell*);
	virtual void inflow_walldist(lexer*,fdm*,ghostcell*,convection*,reini*,ioflow*);
	virtual void discharge(lexer*,fdm*,ghostcell*);
	virtual void inflow(lexer*,fdm*,ghostcell*,field&,field&,field&);
	virtual void rkinflow(lexer*,fdm*,ghostcell*,field&,field&,field&);
	virtual void fsfinflow(lexer*,fdm*,ghostcell*);
	virtual void fsfrkin(lexer*,fdm*,ghostcell*,field&);
	virtual void fsfrkout(lexer*,fdm*,ghostcell*,field&);
	virtual void fsfrkinV(lexer*,fdm*,ghostcell*,vec&);
	virtual void fsfrkoutV(lexer*,fdm*,ghostcell*,vec&);
	virtual void fsfrkinVa(lexer*,fdm*,ghostcell*,vec&);
	virtual void fsfrkoutVa(lexer*,fdm*,ghostcell*,vec&);
    virtual void iogcb_update(lexer*,fdm*,ghostcell*);
    virtual void isource(lexer*,fdm*,ghostcell*,vrans*);
    virtual void jsource(lexer*,fdm*,ghostcell*,vrans*);
    virtual void ksource(lexer*,fdm*,ghostcell*,vrans*);
    virtual void pressure_io(lexer*,fdm*,ghostcell*);
    virtual void turbulence_io(lexer*,fdm*,ghostcell*);
    virtual void veltimesave(lexer*,fdm*,ghostcell*,vrans*);
    virtual void flowfile(lexer*,fdm*,ghostcell*,turbulence*);
    
    virtual void wavegen_precalc(lexer*,ghostcell*);
    virtual void wavegen_precalc_ini(lexer*,ghostcell*);
    virtual void u_relax(lexer*,fdm*,ghostcell*,field&);
    virtual void v_relax(lexer*,fdm*, ghostcell*,field&);
    virtual void w_relax(lexer*,fdm*, ghostcell*,field&);
    virtual void p_relax(lexer*,fdm*,ghostcell*,field&);
	virtual void phi_relax(lexer*,ghostcell*,field&);
    virtual void vof_relax(lexer*,ghostcell*,field&);
    virtual void turb_relax(lexer*,fdm*,ghostcell*,field&);
    virtual void U_relax(lexer*,ghostcell*,double*,double*);
    virtual void V_relax(lexer*,ghostcell*,double*,double*);
    virtual void W_relax(lexer*,ghostcell*,double*,double*);
    virtual void P_relax(lexer*,ghostcell*,double*);
    virtual void WL_relax(lexer*,ghostcell*,slice&,slice&);
    virtual void fi_relax(lexer*,ghostcell*,field&,field&);
    virtual void fivec_relax(lexer*, ghostcell*, double*);
    virtual void fifsf_relax(lexer*, ghostcell*, slice&);
    virtual void visc_relax(lexer*, ghostcell*, slice&);
    virtual void eta_relax(lexer*,ghostcell*,slice&);
    virtual void um_relax(lexer*,ghostcell*,slice&,slice&,slice&);
    virtual void vm_relax(lexer*,ghostcell*,slice&,slice&,slice&);
	virtual void wm_relax(lexer*,ghostcell*,slice&,slice&,slice&);
    virtual void ws_relax(lexer*,ghostcell*,slice&,slice&,slice&);
	virtual void pm_relax(lexer*,ghostcell*,slice&);
    
    virtual void wavegen_2D_precalc(lexer*,fdm2D*,ghostcell*);
    virtual void wavegen_2D_precalc_ini(lexer*,ghostcell*);
    
    virtual void discharge2D(lexer*,fdm2D*,ghostcell*);
    virtual void waterlevel2D(lexer*,fdm2D*,ghostcell*,slice&);
    virtual void Qin2D(lexer*,fdm2D*,ghostcell*);
	virtual void Qout2D(lexer*,fdm2D*,ghostcell*);
    virtual void inflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&);
	virtual void rkinflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&);
	virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);
	virtual void full_initialize2D(lexer*,fdm2D*,ghostcell*);
    
    virtual double wave_fsf(lexer*,ghostcell*,double);
	virtual int iozonecheck(lexer*,fdm*);
    
    virtual void ini(lexer*,fdm*,ghostcell*);
    virtual void ini_fnpf(lexer*,fdm_fnpf*,ghostcell*);
    virtual void inflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,double*,double*,slice&,slice&);
    virtual void rkinflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&);
    virtual void ini2D(lexer*,fdm2D*,ghostcell*);
    virtual void ini_ptf(lexer*,fdm*,ghostcell*);
    
    // nhflow
    virtual void wavegen_precalc_nhflow(lexer*,fdm_nhf*,ghostcell*);
    virtual void wavegen_precalc_ini_nhflow(lexer*,fdm_nhf*,ghostcell*);
    virtual void ini_nhflow(lexer*,fdm_nhf*,ghostcell*);
    virtual void discharge_nhflow(lexer*,fdm_nhf*,ghostcell*);
    virtual void inflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*);
    virtual void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*);
    virtual void isource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*);
    virtual void jsource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*);
    virtual void ksource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*);
    virtual void fsfinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,slice&);
    
    virtual void vrans_sed_update(lexer*,fdm*,ghostcell*,vrans*);
	
	
private:
	double omega_x,omega_y;
	double theta_x,theta_y;
	double dist_x,dist_y,dist_z;
    
    patchBC_interface *pBC;
};

#endif

