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

#ifndef SEDIMENT_F_H_
#define SEDIMENT_F_H_

#include"sediment.h"
#include"sliceint4.h"
#include"slice4.h"
#include"field4a.h"
#include"increment.h"
#include<vector>

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
class nhflow_scalar_convection;
class nhflow_diffusion;
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
class sediment_roughness;
class sediment_exner;
class fdm_nhf;
class fdm2D;
class ioflow;
class sediment_waves;
class bedload_waves;
class sediment_susp2D;
class slice1;
class slice2;
class slice4;

// mesh refinement and multi-block grids: the sediment step of all grids (nhflow_amr, sflow_amr), called by
// the level-0 sediment_f instead of its own step (sediment_amr.cpp)
class sediment_grids
{
public:
    virtual ~sediment_grids(){}
    virtual void sediment_step_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, class sediment_f*) {}
    virtual void sediment_step_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, class sediment_f*) {}
};

using namespace std;

class sediment_f : public sediment, public increment
{
public:
    sediment_f(lexer*,ghostcell*,turbulence*,patchBC_interface*);
    // a grid of the mesh refinement or a block (amrgrid 1): no log, no probes, no messages
    sediment_f(lexer*,ghostcell*,turbulence*,patchBC_interface*,int amrgrid);
	virtual ~sediment_f();

    // checkpoints (P 48, I 45), sediment_f.cpp
    void ckpt_state(lexer*, std::vector<double>&) override;
    bool ckpt_restore(lexer*, ghostcell*, const std::vector<double>&) override;
    void ckpt_rank_state(lexer*, std::vector<double>&) override;
    bool ckpt_rank_restore(lexer*, ghostcell*, const std::vector<std::vector<double>>&, bool) override;
    double* ckpt_wln() override;
    // multi-block grids: the cells beyond the open faces of a block (G 54) carry the sediment flag as the ghost cells
    // of the native in- and outflow (lexer indices (i - imin)*jmax + j - jmin; active_sflow keeps them open)
    std::vector<int> dfbed_open;
    
    // CFD interface
    void start_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*) override;
    void ini_cfd(lexer*,fdm*,ghostcell*) override final;
    void start_susp(lexer*, fdm*, ghostcell*, ioflow*, solver*) override final;
    void update_cfd(lexer*,fdm*,ghostcell*,ioflow*,reinitopo*) override final;
    
    void sediment_logic(lexer*,ghostcell*,turbulence*);
    void sediment_algorithm_cfd(lexer*, fdm*, ghostcell*, ioflow*, reinitopo*, solver*);
    void prep_cfd(lexer*,fdm*,ghostcell*);
    void fill_PQ_cfd(lexer*,fdm*,ghostcell*);
    void active_cfd(lexer*,fdm*,ghostcell*);
    void active_ini_cfd(lexer*,fdm*,ghostcell*);
    void bedchange_update(lexer*, ghostcell*);
    
    // NHFLOW interface
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*) override;
    void ini_nhflow(lexer*, fdm_nhf*, ghostcell*) override final;
    void start_susp_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*, solver*) override final;
    void update_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override final;
    
    void RK2_step1_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override {};
    void RK2_step2_nhflow(lexer*,fdm_nhf*,ghostcell*,ioflow*) override {};
    
    void sediment_algorithm_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    void prep_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void fill_PQ_nhflow(lexer*,fdm_nhf*,ghostcell*);
    void active_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void active_ini_nhflow(lexer*, fdm_nhf*, ghostcell*);
    
    // SFLOW interface
    void start_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&) override;
    void ini_sflow(lexer*, fdm2D*, ghostcell*) override final;
    void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*) override final;
    void sediment_algorithm_sflow(lexer*, fdm2D*, ghostcell*, ioflow*, slice&, slice&);
    void prep_sflow(lexer*, fdm2D*, ghostcell*,slice&,slice&);
    void fill_PQ_sflow(lexer*,fdm2D*,ghostcell*,slice&,slice&);
    void active_sflow(lexer*, fdm2D*, ghostcell*);
    void active_ini_sflow(lexer*, fdm2D*, ghostcell*);
    
    // SFLOW: bed shear stress with the waves (S 110), shared by the sediment step and the suspended load
    void physics_sflow(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    // SFLOW: depth-averaged suspended load (S 120), every flow step once the sediment transport has started
    void susp_sflow(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    // SFLOW wave formulas (S 11 6-9): the face directions of the Exner step from the transport direction
    // (mode 1, s->P and s->Q saved), mode 0 restores the flow velocities
    void transport_direction(lexer*, ghostcell*, int mode);
    void check_sflow_waves(lexer*);
    
    void active_zone(lexer*, ghostcell*);
    void dfbed_comms(lexer*, ghostcell*);
    
    
    // ---

    void ini_parameters(lexer*, ghostcell*);
    void ini_guard(lexer*, ghostcell*);
	
    void relax(lexer*,ghostcell*) override final;
	double bedshear_point(lexer*,ghostcell*) override final;
    
    double qbeval(int,int) override final;
    void qbeget(int,int,double) override final;

    double bedzhval(int,int) override final;

    void ctimesave(lexer*, fdm*) override final;
    
    void fill_bedk(lexer*,fdm*,ghostcell*);
	void bedlevel(lexer*,ghostcell*);
    void waterlevel(lexer*,fdm*,ghostcell*);
	void topo_zh_update(lexer*,fdm*,ghostcell*,sediment_fdm*);
    void volume_calc(lexer*,fdm*,ghostcell*);
    
    // print
    void print_probes(lexer*, ghostcell*,sediment_fdm*, ioflow*) override final;
    void print_particles(lexer*,sediment_fdm*) override final {};
    
    void print_2D_bedload(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_bedload(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_bedload(lexer*,ofstream&) override final;
    void name_ParaView_bedload(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_bedload(lexer*, int*, int &) override final;
    void offset_ParaView_bedload(lexer*, int*, int &) override final;
    
	void print_2D_bedshear(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_bedshear(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_bedshear(lexer*,ofstream&) override final;
    void name_ParaView_bedshear(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_bedshear(lexer*, int*, int &) override final;
    void offset_ParaView_bedshear(lexer*, int*, int &) override final;
    
    void print_2D_parameter1(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_parameter1(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_parameter1(lexer*,ofstream&) override final;
    void name_ParaView_parameter1(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_parameter1(lexer*, int*, int &) override final;
    void offset_ParaView_parameter1(lexer*, int*, int &) override final;
    
    void print_2D_parameter2(lexer*, ghostcell*,ofstream&) override final;
    void print_3D_parameter2(lexer*, ghostcell*, std::vector<char>&, size_t&) override final;
	void name_ParaView_parallel_parameter2(lexer*,ofstream&) override final;
    void name_ParaView_parameter2(lexer*, ostream&, int*, int &) override final;
    void offset_ParaView_2D_parameter2(lexer*, int*, int &) override final;
    void offset_ParaView_parameter2(lexer*, int*, int &) override final;
    
    void print_3D_CPM(lexer*, ghostcell*,  std::vector<char>&, size_t&) override final {};
    void name_ParaView_parallel_CPM(lexer*, ofstream&) override final {};
    void name_ParaView_CPM(lexer*, ostream&, int*, int &) override final {};
    void offset_ParaView_CPM(lexer*, int*, int &) override final {};
    
    
    void log_ini(lexer*);
    void sedimentlog(lexer*);
    
    // mesh refinement and multi-block grids (sediment_amr.cpp): the step of one grid in parts
    //   amr_physics: velocities, bed slope, roughness, bed shear, bed load and its direction, the reference
    //                concentration, the relaxation zones, qb (sediment_exner::amr_qb)
    //   amr_finish:  relaxation zones, the bed to the flow (the sand slide: pslide in steps, reefamr_sed.h)
    void amr_physics_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void amr_physics_sflow(lexer*, fdm2D*, ghostcell*);
    void amr_finish_nhflow(lexer*, fdm_nhf*, ghostcell*, ioflow*);
    void amr_finish_sflow(lexer*, fdm2D*, ghostcell*, ioflow*);
    sediment_exner* exner();
    sediment_grids *grids = nullptr;    // set by nhflow_amr / sflow_amr when they run the sediment step
    int amrgrid = 0;
    sediment_fdm *s = nullptr;
    bedload *pbed = nullptr;  
    bedconc *pcbed = nullptr;
    sandslide *pslide = nullptr;
    topo_relax *prelax = nullptr;
    vrans *pvrans = nullptr;
    bedslope *pslope = nullptr;
    bedshear_reduction *preduce = nullptr;
    topo *ptopo = nullptr;
    suspended *psusp = nullptr;
    nhflow_suspended *pnhfsusp = nullptr;
    nhflow_scalar_convection *pnhfsuspdisc = nullptr;
    nhflow_diffusion *pnhfsuspdiff = nullptr;
    
    diffusion *psuspdiff = nullptr;
    convection *psuspdisc = nullptr;
	bedshear *pbedshear = nullptr;
    patchBC_interface *pBC;
    bedload_direction *pbeddir = nullptr;
    bedprobe_point *pbedpt = nullptr;
	bedprobe_line_x *pbedlinex = nullptr;
	bedprobe_line_y *pbedliney = nullptr;
	bedprobe_max *pbedmax = nullptr;
	bedshear_probe *pbedshearprobe = nullptr;
	bedshear_max *pbedshearmax = nullptr;
    sediment_f *psed;
    ofstream sedlogout;
    sediment_roughness *pks = nullptr; 
    
    // SFLOW: waves (S 110), wave formulas (S 11 6-9), depth-averaged suspended load (S 120)
    sediment_waves *pwaves = nullptr;
    bedload_waves *pbedw = nullptr;
    sediment_susp2D *psusp2D = nullptr;
    slice1 *Psave = nullptr;
    slice2 *Qsave = nullptr;
    slice4 *Vtx = nullptr, *Vty = nullptr;
    
    double starttime;
    
    int volume_token,sedcalc;
    int gcval_eta;
    double volume0;
    
private:
	
};

#endif
