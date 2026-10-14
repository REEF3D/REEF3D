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

#ifndef SEDIMENT_EXNER_H_
#define SEDIMENT_EXNER_H_

#include"topo.h"
#include"slice4.h"
#include"vec2D.h"
#include"matrix2D.h"

class bedconc;
class topo_relax;
class solver2D;
class turbulence;
class ghostcell;

using namespace std;

class sediment_exner final : public topo, public increment
{
public:
	sediment_exner(lexer*, ghostcell*);
	virtual ~sediment_exner();
	void start(lexer*, ghostcell*, sediment_fdm*) override final;
    void start_RK(lexer*, ghostcell*, sediment_fdm*) override final;
    void  timestep(lexer*,ghostcell*,sediment_fdm*) override final;
    void timestep_curv(lexer* p, ghostcell *pgc, sediment_fdm *s);     // curvilinear grid (G 50)
    
    // mesh refinement and multi-block grids (nhflow_amr, sflow_amr; sediment_exner_amr.cpp): the step
    // of start() in parts, so that the grids share the face fluxes of their interfaces and one dtsed
    //   amr_qb:   qb (equilibrium or non-equilibrium), cleared outside the sediment cells, ghost cells
    //   amr_flux: face fluxes per unit length, Fx(i,j) on face i+1/2, Fy(i,j) on face j+1/2 (S 31 1/2,
    //             S 32 1/2; Cartesian and curvilinear grids)
    //   amr_vz:   vz from the fluxes (and the suspended exchange); the filter (S 100) is the driver's
    //             (sediment_xfilter over all grids)
    //   amr_vzmax: the largest |vz| of the cells with mask > 0.5 (nullptr: all), rank-local
    //   amr_bed:  dh = dtsed*vz, the dry deposit, bedzh += dh
    //   amr_dx:   the smallest cell width of the grid (rank-local), for dtsed
    // amr_maxvz and amr_mindx (> 0) take part in timestep(): the maxima and minima of the other grids.
    void amr_qb(lexer*, ghostcell*, sediment_fdm*);
    void amr_flux(lexer*, ghostcell*, sediment_fdm*, slice&, slice&);
    double amr_vz(lexer*, ghostcell*, sediment_fdm*, slice&, slice&, slice*);
    double amr_vzmax(lexer*, sediment_fdm*, slice*);
    void amr_bed(lexer*, ghostcell*, sediment_fdm*);
    double amr_dx(lexer*);
    double amr_maxvz = 0.0, amr_mindx = 0.0;

private:
    void start_mixture(lexer*,ghostcell*,sediment_fdm*);
    
    void topovel1(lexer*,ghostcell*,sediment_fdm*);
    void topovel2(lexer*,ghostcell*,sediment_fdm*);
    void topovel3(lexer*,ghostcell*,sediment_fdm*);
    void topovel_curv(lexer*,ghostcell*,sediment_fdm*,int);   // curvilinear grid (G 50): S 31 1, 2 (sediment_exner_curv.cpp)
    void face_closure(lexer*,double&,double&,double&,double&);
    void qb_clear(lexer*,sediment_fdm*);
    void dry_deposit(lexer*,sediment_fdm*);
    
    double ramp_dt(lexer*);
    void  non_equillibrium_solve(lexer*,ghostcell*,sediment_fdm*);
    void non_equillibrium_solve_curv(lexer* p, ghostcell *pgc, sediment_fdm *s);     // curvilinear grid (G 50)
    double  susp_ED(lexer*,ghostcell*,sediment_fdm*);
    void  susp_qs(lexer*,ghostcell*,sediment_fdm*);
    
    void filter(lexer*,ghostcell*,slice&,int,int);
    void filter_curv(lexer*,ghostcell*,slice&,int,int);     // curvilinear grid (G 50)   
    topo_relax *prelax;
    sediment_exnerdisc *pdx;
    solver2D *psolv;
        vec2D xvec,rhsvec;

	matrix2D M;
    slice4 qbx,qby;
    slice4 qbn;
    slice4 cxn,cyn;
    
	int gcval_topo;
    int noneq_ini;
    int frac_k;    
	double starttime;
    double maxdh,maxvz,maxdhnet,maxvznet;
    double dtsed0;
	double vz;
	double vzmax;
    double rhosed, rhowat, g, d50;
    double Ls;
    double tau_eff, shearvel_eff, shields_eff;
    double tau_crit, shearvel_crit, shields_crit;
    
    slice4 q0;
    
    double maxvz_k=0.0;   // multi-fraction bed: max |vz| of the single fractions (time step)
    double simtime_last=-1.0;   // flow time of the previous bed update
    double morfac=0.0, c_U=0.0; // effective morphological factor S35*dtsed/flow time, bed celerity/flow speed
};

#endif

