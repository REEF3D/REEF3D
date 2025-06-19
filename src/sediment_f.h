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

#ifndef SEDIMENT_F_H_
#define SEDIMENT_F_H_

#include"sediment.h"
#include"sliceint4.h"
#include"slice4.h"
#include"field4a.h"
#include"increment.h"

class bedload;
class bedconc;
class sandslide;
class topo_relax;
class bedshear;
class vrans;
class turbulence;
class sediment_fdm;
class bedslope;
class bedshear_reduction;
class suspended;
class nhflow_suspended;
class diffusion;
class convection;
class patchBC_interface;
class bedload_direction;
class bedprobe_point;
class bedprobe_max;
class bedshear_probe;
class bedshear_max;
class bedprobe_line_x;
class bedprobe_line_y;

using namespace std;

class sediment_f : public sediment, public increment
{
public:
    sediment_f(lexer*,fdm*,ghostcell*,turbulence*, patchBC_interface*);
	virtual ~sediment_f();
    
    // CFD interface
    virtual void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*);
    virtual void ini_cfd(lexer*,fdm*,ghostcell*);
    virtual void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*);
    virtual void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*);
    
    void sediment_logic(lexer*,fdm*,ghostcell*,turbulence*);
    void sediment_algorithm_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*);
    void prep_cfd(lexer*,fdm*,ghostcell*);
    void fill_PQ_cfd(lexer*,fdm*,ghostcell*);
    void active_cfd(lexer*,fdm*,ghostcell*);
    void active_ini_cfd(lexer*,fdm*,ghostcell*);
    void bedchange_update(lexer*, ghostcell*);
    
    // NHFLOW interface
    virtual void start_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    virtual void ini_nhflow(lexer*, fdm_nhf*, ghostcell*);
    virtual void start_susp_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*);
    virtual void update_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*);
    
    void sediment_algorithm_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    void prep_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void fill_PQ_nhflow(lexer*,fdm_nhf*,ghostcell*);
    void active_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void active_ini_nhflow(lexer*, fdm_nhf*, ghostcell*);
    
    // SFLOW interface
    virtual void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&);
    virtual void ini_sflow(lexer*, fdm2D*, ghostcell*);
    virtual void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*);
    void sediment_algorithm_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&);
    void prep_sflow(lexer*, fdm2D*, ghostcell*,slice&,slice&);
    void fill_PQ_sflow(lexer*,fdm2D*,ghostcell*,slice&,slice&);
    void active_sflow(lexer*, fdm2D*, ghostcell*);
    void active_ini_sflow(lexer*, fdm2D*, ghostcell*);
    
    
    // ---

    virtual void ini_parameters(lexer*, ghostcell*);
    virtual void ini_guard(lexer*, ghostcell*);
	
    virtual void relax(lexer*,ghostcell*);
	virtual double bedshear_point(lexer*,ghostcell*);
    
    virtual double qbeval(int,int);
    virtual void qbeget(int,int,double);

    virtual double bedzhval(int,int);

    virtual void ctimesave(lexer*, fdm*);
    
    void fill_bedk(lexer*,fdm*,ghostcell*);
	void bedlevel(lexer*,ghostcell*);
    void waterlevel(lexer*,fdm*,ghostcell*);
	void topo_zh_update(lexer*,fdm*,ghostcell*,sediment_fdm*);
    void volume_calc(lexer*,fdm*,ghostcell*);
	void filter(lexer*,ghostcell*,slice&,int,int);
    
    // print
    virtual void print_probes(lexer*, ghostcell*,sediment_fdm*, ioflow*);
    
    virtual void print_2D_bedload(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_bedload(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedload(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_bedload(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &);
    
	virtual void print_2D_bedshear(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_bedshear(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedshear(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_2D_parameter1(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_parameter1(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameter1(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_parameter1(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_2D_parameter2(lexer*, ghostcell*,ofstream&);
    virtual void print_3D_parameter2(lexer*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameter2(lexer*, ghostcell*,ofstream&);
    virtual void name_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtp_parameter2(lexer*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &);
    

private:
    
    void log_ini(lexer*);
    void sedimentlog(lexer*);
    sediment_fdm *s;
    bedload *pbed;  
    bedconc *pcbed;
    sandslide *pslide;
    topo_relax *prelax;
    vrans *pvrans;
    bedslope *pslope;
    bedshear_reduction *preduce;
    topo *ptopo;
    suspended *psusp;
    nhflow_suspended *pnhfsusp;
    diffusion *psuspdiff;
    convection *psuspdisc;
	bedshear *pbedshear;
    patchBC_interface *pBC;
    bedload_direction *pbeddir;
    bedprobe_point *pbedpt;
	bedprobe_line_x *pbedlinex;
	bedprobe_line_y *pbedliney;
	bedprobe_max *pbedmax;
	bedshear_probe *pbedshearprobe;
	bedshear_max *pbedshearmax;
    sediment_f *psed;
    ofstream sedlogout;
    
    double starttime;
    
    int volume_token,sedcalc;
    int gcval_eta;
    double volume0;
	
};

#endif
