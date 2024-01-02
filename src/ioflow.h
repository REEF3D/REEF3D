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

class lexer;
class fdm;
class ghostcell;
class field;
class vec;
class convection;
class reini;
class vrans;
class slice;
class fdm2D;
class fdm_fnpf;
class fdm_nhf;
class turbulence;
class patchBC_interface;

using namespace std;

#ifndef IOFLOW_H_
#define IOFLOW_H_

class ioflow
{
public:
    virtual void gcio_update(lexer*,fdm*,ghostcell*)=0;
	virtual void inflow_walldist(lexer*,fdm*,ghostcell*,convection*,reini*,ioflow*)=0;
	virtual void discharge(lexer*,fdm*,ghostcell*)=0;
	virtual void inflow(lexer*,fdm*,ghostcell*,field&,field&,field&)=0;
	virtual void rkinflow(lexer*,fdm*,ghostcell*,field&,field&,field&)=0;
	virtual void fsfinflow(lexer*,fdm*,ghostcell*)=0;
	virtual void fsfrkin(lexer*,fdm*,ghostcell*,field&)=0;
	virtual void fsfrkout(lexer*,fdm*,ghostcell*,field&)=0;
	virtual void fsfrkinV(lexer*,fdm*,ghostcell*,vec&)=0;
	virtual void fsfrkoutV(lexer*,fdm*,ghostcell*,vec&)=0;
	virtual void fsfrkinVa(lexer*,fdm*,ghostcell*,vec&)=0;
	virtual void fsfrkoutVa(lexer*,fdm*,ghostcell*,vec&)=0;
    virtual void iogcb_update(lexer*,fdm*,ghostcell*)=0;
    virtual void isource(lexer*,fdm*,ghostcell*,vrans*)=0;
    virtual void jsource(lexer*,fdm*,ghostcell*,vrans*)=0;
    virtual void ksource(lexer*,fdm*,ghostcell*,vrans*)=0;
    virtual void pressure_io(lexer*,fdm*,ghostcell*)=0;
    virtual void turbulence_io(lexer*,fdm*,ghostcell*)=0;
    virtual void veltimesave(lexer*,fdm*,ghostcell*,vrans*)=0;
    
    virtual void flowfile(lexer*,fdm*,ghostcell*,turbulence*)=0;
    
    
    virtual void wavegen_precalc(lexer*,ghostcell*)=0;
    virtual void wavegen_precalc_ini(lexer*,ghostcell*)=0;
    virtual void u_relax(lexer*,fdm*,ghostcell*,field&)=0;
    virtual void v_relax(lexer*,fdm*,ghostcell*,field&)=0;
    virtual void w_relax(lexer*,fdm*,ghostcell*,field&)=0;
    virtual void p_relax(lexer*,fdm*,ghostcell*,field&)=0;
	virtual void phi_relax(lexer*,ghostcell*,field&)=0;
    virtual void vof_relax(lexer*,ghostcell*,field&)=0;
    virtual void turb_relax(lexer*,fdm*,ghostcell*,field&)=0;
    virtual void U_relax(lexer*,ghostcell*,double*,double*)=0;
    virtual void V_relax(lexer*,ghostcell*,double*,double*)=0;
    virtual void W_relax(lexer*,ghostcell*,double*,double*)=0;
    virtual void P_relax(lexer*,ghostcell*,double*)=0;
    virtual void WL_relax(lexer*,ghostcell*,slice&,slice&)=0;
    virtual void fi_relax(lexer*,ghostcell*,field&,field&)=0;
    virtual void fivec_relax(lexer*, ghostcell*, double*)=0;
    virtual void fifsf_relax(lexer*, ghostcell*, slice&)=0;
    virtual void visc_relax(lexer*, ghostcell*, slice&)=0;
    virtual void eta_relax(lexer*,ghostcell*,slice&)=0;
    virtual void um_relax(lexer*,ghostcell*,slice&,slice&,slice&)=0;
    virtual void vm_relax(lexer*,ghostcell*,slice&,slice&,slice&)=0;
	virtual void wm_relax(lexer*,ghostcell*,slice&,slice&,slice&)=0;
    virtual void ws_relax(lexer*,ghostcell*,slice&,slice&,slice&)=0;
	virtual void pm_relax(lexer*,ghostcell*,slice&)=0;
    
    virtual void wavegen_2D_precalc(lexer*,fdm2D*,ghostcell*)=0;
    virtual void wavegen_2D_precalc_ini(lexer*,ghostcell*)=0;
    
    virtual void discharge2D(lexer*,fdm2D*,ghostcell*)=0;
    virtual void waterlevel2D(lexer*,fdm2D*,ghostcell*,slice&)=0;
    virtual void Qin2D(lexer*,fdm2D*,ghostcell*)=0;
	virtual void Qout2D(lexer*,fdm2D*,ghostcell*)=0;
    virtual void inflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&)=0;
	virtual void rkinflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&)=0;
	virtual void isource2D(lexer*,fdm2D*,ghostcell*)=0;
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*)=0;
	virtual void full_initialize2D(lexer*,fdm2D*,ghostcell*)=0;
    
    virtual void ini(lexer*,fdm*,ghostcell*)=0;
    
    virtual void ini_fnpf(lexer*,fdm_fnpf*,ghostcell*)=0;
    virtual void inflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,double*,double*,slice&,slice&)=0;
    virtual void rkinflow_fnpf(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&)=0;
    
    virtual void ini_ptf(lexer*,fdm*,ghostcell*)=0;
    
    // nhflow
    virtual void wavegen_precalc_nhflow(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void wavegen_precalc_ini_nhflow(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void ini_nhflow(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void discharge_nhflow(lexer*,fdm_nhf*,ghostcell*)=0;
    virtual void inflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*)=0;
    virtual void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*)=0;
    virtual void isource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*)=0;
    virtual void jsource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*)=0;
    virtual void ksource_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*)=0;
    virtual void fsfinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,slice&)=0;


    virtual void ini2D(lexer*,fdm2D*,ghostcell*)=0;

    virtual double wave_fsf(lexer*,ghostcell*,double)=0;
	
	virtual int iozonecheck(lexer*,fdm*)=0;
    
    virtual void vrans_sed_update(lexer*,fdm*,ghostcell*,vrans*)=0;
	
};

#endif
