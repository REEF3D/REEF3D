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

#ifndef SEDIMENT_PART2_H_
#define SEDIMENT_PART2_H_

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
class bedshear_reduction;
class suspended;
class diffusion;
class convection;
class patchBC_interface;
class bedload_direction;
class bedslope;
class partres;
using namespace std;

class sediment_part : public sediment, public increment
{
public:
    sediment_part(lexer*,fdm*,ghostcell*,turbulence*, patchBC_interface*);
	virtual ~sediment_part();
    
    // CFD interface
    void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*) override;
    void ini_cfd(lexer*,fdm*,ghostcell*) override;
    void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*) override {};
    void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*) override;
    
    void sediment_logic(lexer*,fdm*,ghostcell*,turbulence*);
    void sediment_algorithm_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, turbulence*);
    void prep_cfd(lexer*,fdm*,ghostcell*){};
    void fill_PQ_cfd(lexer*,fdm*,ghostcell*){};
    void active_cfd(lexer*,fdm*,ghostcell*){};
    void active_ini_cfd(lexer*,fdm*,ghostcell*){};
    void bedchange_update(lexer*, ghostcell*){};
    
    // NHFLOW interface
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*) override {};
    void ini_nhflow(lexer*, fdm_nhf*, ghostcell*) override {};
    void start_susp_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*) override {};
    void update_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override {};
    
    void sediment_algorithm_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*){};
    void prep_nhflow(lexer*, fdm_nhf*, ghostcell*){};
    void fill_PQ_nhflow(lexer*,fdm_nhf*,ghostcell*){};
    void active_nhflow(lexer*, fdm_nhf*, ghostcell*){};
    void active_ini_nhflow(lexer*, fdm_nhf*, ghostcell*){};
    
    // SFLOW interface
    void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&) override {};
    void ini_sflow(lexer*, fdm2D*, ghostcell*) override {};
    void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*) override {};
    void sediment_algorithm_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&){};
    void prep_sflow(lexer*, fdm2D*, ghostcell*,slice&,slice&){};
    void fill_PQ_sflow(lexer*,fdm2D*,ghostcell*,slice&,slice&){};
    void active_sflow(lexer*, fdm2D*, ghostcell*){};
    void active_ini_sflow(lexer*, fdm2D*, ghostcell*){};
    
    
    // ---


    void relax(lexer*,ghostcell*) override {};
	double bedshear_point(lexer*,ghostcell*) override {};
    
    double qbeval(int,int) override {};
    void qbeget(int,int,double) override {};

    double bedzhval(int,int) override {};

    void ctimesave(lexer*, fdm*) override {};
    
    void fill_bedk(lexer*,fdm*,ghostcell*){};
	void bedlevel(lexer*,fdm*,ghostcell*){};
    void waterlevel(lexer*,fdm*,ghostcell*){};
	void topo_zh_update(lexer*,fdm*,ghostcell*,sediment_fdm*){};
    void volume_calc(lexer*,fdm*,ghostcell*){};
	void filter(lexer*,ghostcell*,slice&,int,int){};
    
    // print
    void print_probes(lexer*, ghostcell*,sediment_fdm*, ioflow*) override {};
    
    void print_2D_bedload(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_bedload(lexer*, ghostcell*,ofstream&) override {};
	void name_pvtu_bedload(lexer*, ghostcell*,ofstream&) override {};
    void name_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtp_bedload(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtu_bedload(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    
	void print_2D_bedshear(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_bedshear(lexer*, ghostcell*,ofstream&) override {};
	void name_pvtu_bedshear(lexer*, ghostcell*,ofstream&) override {};
    void name_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtp_bedshear(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtu_bedshear(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    
    void print_2D_parameter1(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_parameter1(lexer*, ghostcell*,ofstream&) override {};
	void name_pvtu_parameter1(lexer*, ghostcell*,ofstream&) override {};
    void name_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtp_parameter1(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtu_parameter1(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    
    void print_2D_parameter2(lexer*, ghostcell*,ofstream&) override {};
    void print_3D_parameter2(lexer*, ghostcell*,ofstream&) override {};
	void name_pvtu_parameter2(lexer*, ghostcell*,ofstream&) override {};
    void name_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtp_parameter2(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    void offset_vtu_parameter2(lexer*, ghostcell*,ofstream&, int*, int &) override {};
    
    
    partres *pst;
    
    field4a por, d50;

private:
    
    void log_ini(lexer*);
    void sedimentlog(lexer*);
    sediment_fdm *s;
    bedload *pbed;  
    vrans *pvrans;
    bedconc *pcbed;
    sandslide *pslide;
    topo_relax *prelax;
    bedshear_reduction *preduce;
    topo *ptopo;
    suspended *psusp;
    diffusion *psuspdiff;
    convection *psuspdisc;
	bedshear *pbedshear;
    patchBC_interface *pBC;
    bedload_direction *pbeddir;
    bedslope *pslope;
    turbulence *pturb;

    ofstream sedlogout;
    
    double starttime;
    
    int volume_token,sedcalc;
    int gcval_eta;
    double volume0;
	
};

#endif
