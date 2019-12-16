/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
#include"resize.h"
#include"increment.h"
#include"field4.h"

class vrans;

using namespace std;

#ifndef IOFLOW_F_H_
#define IOFLOW_F_H_


class ioflow_f : public ioflow, private resize_class, public increment
{

public:
	ioflow_f(lexer*, ghostcell*);
	virtual ~ioflow_f();
	virtual void gcio_update(lexer*,fdm*,ghostcell*);
	virtual void inflow_walldist(lexer*,fdm*,ghostcell*,convection*,reini*,ioflow*);
	virtual void discharge(lexer*,fdm*,ghostcell*);
	virtual void Qin(lexer*,fdm*,ghostcell*);
	virtual void Qout(lexer*,fdm*,ghostcell*);
	virtual void inflow(lexer*,fdm*,ghostcell*,field&,field&,field&);
	virtual void rkinflow(lexer*,fdm*,ghostcell*,field&,field&,field&);
	virtual void inflow_plain(lexer*,fdm*,ghostcell*,field&,field&,field&);
    virtual void inflow_log(lexer*,fdm*,ghostcell*,field&,field&,field&);
    virtual void outflow_log(lexer*,fdm*,ghostcell*,field&,field&,field&);
    virtual void outflow_plain(lexer*,fdm*,ghostcell*,field&,field&,field&);
    virtual void outflow_water(lexer*,fdm*,ghostcell*,field&,field&,field&);
    virtual void inflow_water(lexer*,fdm*,ghostcell*,field&,field&,field&);
	virtual void fsfinflow(lexer*,fdm*,ghostcell*);
	virtual void fsfrkin(lexer*,fdm*,ghostcell*,field&);
	virtual void fsfrkout(lexer*,fdm*,ghostcell*,field&);
	virtual void fsfrkinV(lexer*,fdm*,ghostcell*,vec&);
	virtual void fsfrkoutV(lexer*,fdm*,ghostcell*,vec&);
	virtual void fsfrkinVa(lexer*,fdm*,ghostcell*,vec&);
	virtual void fsfrkoutVa(lexer*,fdm*,ghostcell*,vec&);
	virtual void iogcb_update(lexer*,fdm*,ghostcell*);
	virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    virtual void pressure_io(lexer*,fdm*,ghostcell*);
    virtual void turbulence_io(lexer*,fdm*,ghostcell*);
    virtual void veltimesave(lexer*,fdm*,ghostcell*);
    virtual void flowfile(lexer*,fdm*,ghostcell*,turbulence*);
    
    virtual void wavegen_precalc(lexer*,ghostcell*);
    virtual void wavegen_precalc_ini(lexer*,ghostcell*);
    virtual void u_relax(lexer*,fdm*,ghostcell*,field&);
    virtual void v_relax(lexer*,fdm*, ghostcell*,field&);
    virtual void w_relax(lexer*,fdm*, ghostcell*,field&);
    virtual void p_relax(lexer*,fdm*,ghostcell*,field&);
	virtual void phi_relax(lexer*,ghostcell*,field&);
    virtual void vof_relax(lexer*,ghostcell*,field&);
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
    virtual void inflow_fnpf(lexer*,ghostcell*,double*,double*,slice&,slice&);
    virtual void ini2D(lexer*,fdm2D*,ghostcell*);
    
    virtual void ini_nhflow(lexer*,fdm*,ghostcell*);
    virtual void nhflow_inflow(lexer*,fdm*,ghostcell*,field&,field&,field&);
    
    virtual void vrans_sed_update(lexer*,fdm*,ghostcell*);
	
	void hydrograph_in_read(lexer*,ghostcell*);
	void hydrograph_out_read(lexer*,ghostcell*);
	double hydrograph_ipol(lexer*,ghostcell*,double**,int);

    void pressure_outlet(lexer*,fdm*,ghostcell*);
    void pressure_inlet(lexer*,fdm*,ghostcell*);
    void pressure_wall(lexer*,fdm*,ghostcell*);
    void pressure_bed(lexer*,fdm*,ghostcell*);
    double local_fsf(lexer*,fdm*,ghostcell*);

    void fsfdistance(lexer*,fdm*,ghostcell*);
    double r1(lexer*,double);
	double r3(double,double);

    void velini(lexer*,fdm*,ghostcell*);

private:

    int n,count;
    double area,Ai,Ao,Ui,fac;
    double zval;
    double *walldin, *walldout;
	int walldin_size, walldout_size;
	
	double distcalc(lexer*,double, double, double);
	double r1(lexer*, double, double);

	double *tan_betaB67,*tan_betaB68,*tan_betaB69,*tan_betaB70,*tan_betaB71;
	double *betaB67,*betaB68,*betaB69,*betaB70,*betaB71;
	double *dist_B67, *dist_B68, *dist_B69,*dist_B70,*dist_B71;
	
	double **hydro_in,**hydro_out;
    int hydro_in_count,hydro_out_count;
	
	double kinval, epsval, eddyval, val;
    
    double Apor,Bpor,porval,partval;
	
	vrans *pvrans;

};

#endif
