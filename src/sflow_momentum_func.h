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

#ifndef SFLOW_MOMENTUM_FUNC_H_
#define SFLOW_MOMENTUM_FUNC_H_

#include"sflow_momentum.h"
#include"increment.h"
#include"slice4.h"

class slice;
class ioflow;
class solver2D;
class sixdof;
class sflow_HLL;
class sflow_fsf;
class sflow_signal_speed;
class sflow_reconstruct;
class sflow_diffusion;
class sflow_pressure;
class sflow_forcing;
class sflow_roughness;
class sflow_rheology;

using namespace std;

// shared functions of the SFLOW HLL time integrators
// prognostic: WL, UH, VH, WH (cell centres)
// diagnostic: U, V, W (cell centres), P, Q (face velocities), ws, hx, hy, hp
//
// SSP Runge-Kutta stage:  q_out = a*q_n + (1-a)*(q_s + dt*L(q_s)),  alpha = 1-a

class sflow_momentum_func : public sflow_momentum, public increment
{
public:
	sflow_momentum_func(lexer*, fdm2D*, ghostcell*, sflow_HLL*, sflow_signal_speed*, sflow_reconstruct*, 
                        sflow_diffusion*, sflow_pressure*, solver2D*, solver2D*, ioflow*, sflow_fsf*, sflow_forcing*, sixdof*);
	virtual ~sflow_momentum_func();
    
    void ini(lexer*, fdm2D*, ghostcell*) override final;
    void stage(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, slice&, slice&, slice&, slice&, double, int, bool);
    
    void reconstruct(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&);
    void velcalc(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, int);
    void ghostcells(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&);
    void face_velocities(lexer*, fdm2D*, ghostcell*, slice&);
    void inflow(lexer*, fdm2D*, ghostcell*, ioflow*);
    
	void irhs(lexer*, fdm2D*);
	void jrhs(lexer*, fdm2D*);
	void krhs(lexer*, fdm2D*);
    void clearrhs(lexer*, fdm2D*);
    
	int gcval_u, gcval_v, gcval_w;
    int gcval_uh, gcval_vh, gcval_wh;
    int gcval_eta;
    int inflow_flag, outflow_flag;
    double starttime;
    
protected:
    sflow_HLL *phll;
    sflow_signal_speed *pss;
    sflow_reconstruct *precon;
	sflow_diffusion *pdiff;
	sflow_pressure *ppress;
	solver2D *psolv;
    solver2D *ppoissonsolv;
	ioflow *pflow;
	sflow_fsf *pfsf;
    sflow_forcing *psfdf;
    sixdof *p6dof;
    sflow_roughness *prough;
    sflow_rheology *prheo;
    
    slice4 UHDIFF,VHDIFF,WHDIFF;
    slice4 Un,Vn,etaS;
    
private:
    int q;
    
    void mpi4(lexer*, ghostcell*, slice&);
    
    void mirror(lexer*, slice&, int, double);
    void neumann(lexer*, slice&, int);
};

#endif
