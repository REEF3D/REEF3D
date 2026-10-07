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

#ifndef NHFLOW_MOMENTUM_FUNC_H_
#define NHFLOW_MOMENTUM_FUNC_H_

#include"nhflow_momentum.h"
#include"nhflow_bcmom.h"
#include"nhflow_sigma.h"

class nhflow_fsf;
class nhflow_signal_speed;
class nhflow_reconstruct;
class nhflow_fsf_reconstruct;
class nhflow_momentum_func;
class nhflow_forcing;
class sixdof;

using namespace std;

// the objects one RK stage works with (the arguments of nhflow_momentum::start)
struct nhflow_stage_obj
{
    ioflow *pflow;
    nhflow_signal_speed *pss;
    nhflow_reconstruct *precon;
    nhflow_convection *pconvec;
    nhflow_diffusion *pdiff;
    nhflow_pressure *ppress;
    solver *ppoissonsolv;
    solver *psolv;
    nhflow *pnhf;
    nhflow_fsf *pfsf;
    nhflow_turbulence *pturb;
    vrans_nhflow *pvrans;
};

// runs the time step of all grids (mesh refinement, nhflow_amr) instead of start()
class nhflow_stage_runner
{
public:
    virtual void step(lexer*, fdm_nhf*, ghostcell*, nhflow_momentum_func*, nhflow_stage_obj&) = 0;
};

class nhflow_momentum_func : public nhflow_momentum, public nhflow_bcmom, public nhflow_sigma
{
public:
	nhflow_momentum_func(lexer*, fdm_nhf*, ghostcell*);
	virtual ~nhflow_momentum_func();
    

    void inidisc(lexer*, fdm_nhf*, ghostcell*, nhflow_fsf*) override final;
    void reconstruct(lexer*, fdm_nhf*, ghostcell*, nhflow_fsf*, nhflow_signal_speed*, nhflow_reconstruct*,slice&,double*,double*,double*,double*,double*,double*);
    void velcalc(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,slice&,double);
    
	void irhs(lexer*,fdm_nhf*,ghostcell*);
	void jrhs(lexer*,fdm_nhf*,ghostcell*);
	void krhs(lexer*,fdm_nhf*,ghostcell*);
    
    void clearrhs(lexer*,fdm_nhf*,ghostcell*);
    
    // one time step as phases, so that the AMR module can interleave the grids (nhflow_amr):
    //   step_begin  in- and outflow of the stage arrays
    //   phase_F     sigma, reconstruction, continuity flux, water level of stage s, omega
    //   phase_M     momentum fluxes and the RK update of UH, VH, WH
    //   phase_P     velocities, forcing, pressure projection (phase_P1, ppress->start, phase_P2)
    //   phase_E     relaxation zones, ghost cells (RK2: sediment, depth update)
    // start() runs them in this order; without mesh refinement the operations are unchanged
    virtual int stages() const = 0;
    virtual void step_begin(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&) = 0;
    virtual void phase_F(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) = 0;
    virtual void phase_M(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) = 0;
    // the implicit diffusion of stage s alone, component m (0 UH, 1 VH, 2 WH), with the arguments
    // phase_M gives it: nhflow_amr (G 31 1) assembles it on all grids for one composite solve and
    // then runs phase_M with a diffusion object that keeps the result
    virtual void phase_D(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int,int) = 0;
    void phase_P(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int);
    virtual void phase_P1(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) = 0;   // velocities, forcing
    virtual void phase_P2(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) = 0;   // velocities, reforcing
    virtual void phase_E(lexer*,fdm_nhf*,ghostcell*,nhflow_stage_obj&,int) = 0;
    // stage values of stage s: water level and UH, VH, WH (s = stages()-1: d->WL, d->UH, ...)
    virtual slice& stage_WL(fdm_nhf*,int) = 0;
    virtual double* stage_UH(fdm_nhf*,int,int) = 0;
    virtual double stage_alpha(int) const = 0;
    
    void attach_runner(nhflow_stage_runner *a) override { prun = a; }
    nhflow_stage_runner *prun = nullptr;
    
    // membranes (X 330, membrane.dat 'coupling iterated'): phase_P repeats the projection (nhflow_forcing::projection)
    nhflow_forcing *pmfrc = nullptr;
    sixdof *pm6dof = nullptr;
	
    

	int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
    
    double starttime;

};

#endif
