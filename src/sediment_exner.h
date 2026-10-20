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
    // the cell fields of the Exner step that the cells around the grids take from their sources: qb, and
    // the bend fields qb/f(theta) (S 86), tan(delta) (S 95), the bed of the slope flux and its open flag
    // (P7) and the lagged fields the upwind relaxation of a grid starts from at its edges: the spiral-flow
    // intensity (S 95 2) and the suspended load (S 64)
    static const int amr_nfield = 7;
    double& amr_field(sediment_fdm*, int n, int ii, int jj);
    slice4& bend_slice(int n) { return (n==1) ? bK : (n==2 ? btd : (n==3 ? bz : (n==4 ? bg : (n==5 ? Isp : qsn)))); }
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
    
    // bend morphology (sediment_exner_bend.cpp): the transverse bed slope effect (S 86) as a face flux
    // -qb/f(theta) dz/dn, and the spiral-flow correction (S 95) of the transport direction, tan(delta)
    // = alpha_I h/R_s (turned towards the centre of the streamline curvature)
    void bend_prep(lexer*,ghostcell*,sediment_fdm*);
    double bend_flux(lexer*,sediment_fdm*,int,int,int);
    double bend_div(lexer*,sediment_fdm*,int,int);
    bool bend_open(lexer*,int,int);
    double bend_td(lexer*,int,int,int);
    static double bend_turn(double un, double ut, double mag, double td);
    slice4 bK,btd,bu,bv,bz,bg;
    // lags of the bend morphology (sediment_exner_lag.cpp): the steady relaxation f + L s.grad(f) = f_e along the
    // flow, the spiral-flow intensity (S 95 2: Ie its equilibrium, Isp the lagged one), the suspended load (S 64:
    // qse its capacity, qsn the lagged one, qbb the bed load), the lengths lss and lsa
    int relax_along(lexer*,ghostcell*,sediment_fdm*,slice4&,slice4&,slice4&,int&);
public:
    void susp_split(lexer*,ghostcell*,sediment_fdm*);
    void lag_refresh(lexer*,ghostcell*,sediment_fdm*,double*,double*);
    // the tolerance of the relaxations and of the passes across the grids (relative to the largest equilibrium)
    static constexpr double lag_tol = 1.0e-6;
    // mesh refinement and multi-block grids: the cells the relaxation holds (> 0: the cells around the grid, which
    // take their values from their sources, and the covered cells, which take the mean of their children)
    slice4 lagfix;
private:
    slice4 lagval;
    void lag_hold(lexer*,ghostcell*,slice4&,int);
    void susp_apply(lexer*,ghostcell*,sediment_fdm*);
    void spiral_lag(lexer*,ghostcell*,sediment_fdm*);
    void spiral_apply(lexer*,ghostcell*,sediment_fdm*);
    void bend_K(lexer*,ghostcell*,sediment_fdm*);
    slice4 Ie,Isp,qse,qsn,qbb,lss,lsa;
    int spiral_ini=0, susp_ini=0, lag_print=1;
    double lag_d0=0.0, lag_f0=0.0, lag_d[2]={0.0,0.0}, lag_f[2]={0.0,0.0};
    // the largest qb/f (1/dx^2 + 1/dy^2) of the grid (rank-local): the slope flux is an explicit diffusion of the
    // bed, timestep() keeps dtsed below its stability limit; amr_bendrate: that of the other grids
public:
    double bend_rate = 0.0, amr_bendrate = 0.0;
private:
    int p_S95=0;

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

