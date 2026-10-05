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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#ifndef SIXDOF_OBJ_NHFLOW_H_
#define SIXDOF_OBJ_NHFLOW_H_

#include"6DOF_obj_2D.h"

class nhflow_reinidisc_fsf;

//  Fluid access of the load models (ship module) in NHFLOW: velocity at points, MPI-summed
class sixdof_fluid_nhflow : public sixdof_fluid
{
public:
    sixdof_fluid_nhflow(lexer *pp, fdm_nhf *dd, ghostcell *gc) : p(pp), d(dd), pgc(gc) {}
    void velocity(int, const double*, double*) override;
private:
    lexer *p;
    fdm_nhf *d;
    ghostcell *pgc;
};

//  6DOF body coupled to REEF3D::NHFLOW: level set of the body on the sigma grid, direct
//  forcing, hull pressure and viscous loads, membranes (X 330), porous bodies (X 16) and
//  the wavemaker objects (X 10 4). The ship-wave mode (X 10 3) uses the 2D pressure patch.

class sixdof_obj_nhflow : public sixdof_obj_2D
{
public:
    
    sixdof_obj_nhflow(lexer*, ghostcell*, int);
	virtual ~sixdof_obj_nhflow();
    
    void initialize_nhflow(lexer*,fdm_nhf*,ghostcell*);
    void initialize_wavemaker(lexer*,fdm_nhf*,ghostcell*,slice&,slice&);
    void hydrodynamic_forces_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&,bool);
    void update_position_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&, bool);
    void update_wavemaker_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&, bool);
    void solve_eqmotion_nhflow(lexer*,fdm_nhf*,ghostcell*,int,bool);
    void solve_eqmotion_oneway_nhflow(lexer*,ghostcell*,int,bool);
    void update_forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, slice&, int);
    void update_forcing_nhflow_wavemaker(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, slice&, int);
    double Hsolidface_nhflow(lexer*, fdm_nhf*, int,int,int);
    double Hpsi_nhflow(lexer*, fdm_nhf*);
    void membrane_forcing_nhflow(lexer*,fdm_nhf*,ghostcell*,double,double*,double*,double*,slice&);
    void membrane_reaction_nhflow(lexer*,fdm_nhf*,ghostcell*,double,slice&,bool);
    bool membrane_iterated();
    void membrane_reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,double,double*,double*,double*,slice&);
    bool membrane_couple_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double,slice&,int);
    void membrane_stabilisation(lexer*,int);
    Eigen::Vector3d umem_n_=Eigen::Vector3d::Zero(), amem_n_=Eigen::Vector3d::Zero();
    double tmem_n_=-1.0;
    Eigen::Vector3d wmem_n_=Eigen::Vector3d::Zero(), almem_n_=Eigen::Vector3d::Zero();
    void update_forcing_nhflow_porous(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, double*, double*, double*, slice&, int);
    void porosity_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void porous_damping_nhflow(lexer*, int);
    void ray_cast_nhflow_grid(lexer*, fdm_nhf*, ghostcell*, int*, int*, int*, double);
    
    // actuator disks of the load models (body-force propellers): momentum source of the
    // component comp (0 x, 1 y, 2 z) added to the right-hand side F (NHFLOW form, times WL)
    void actuator_source(lexer*, fdm_nhf*, ghostcell*, slice&, int, double*);
    bool actuator_warned = false;
    double nhflow_dsm() const {return DSM;}
    struct nhflow_grid { lexer *p; fdm_nhf *d; slice *WL; };
    std::function<nhflow_grid(double,double)> amr_grid_nhflow;
    // mesh refinement with placed patches (G 40 1): this rank takes the triangle with the centroid
    // (x,y) (it holds the finest grid there), instead of the level-0 subdomain test
    std::function<bool(double,double)> amr_owner_nhflow;

private:
    
    void hydrodynamic_forces_nhflow_volume(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*, slice&, int, bool);
    void externalForces_nhflow(lexer*, fdm_nhf*, ghostcell*, double, bool);
    void netForces_nhflow(lexer*, fdm_nhf*, ghostcell*, double, bool);
    void geometry_parameters_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void geometry_ls_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void update_trimesh_nhflow(lexer*, fdm_nhf*, ghostcell*, bool);
    void ray_cast(lexer*, fdm_nhf*, ghostcell*);
    int  clip_facet_poly(lexer*,double,double,double,double,double,double,double,double,double, double,double*,double*,double*);
    void nhflow_reini_RK2(lexer*, fdm_nhf*, ghostcell*, double*);
    void allocate(lexer*,fdm_nhf*,ghostcell*);
    void deallocate(lexer*,fdm_nhf*,ghostcell*);
    void print_force(lexer*,fdm_nhf*,ghostcell*);
    void forces_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void force_calc_stl(lexer*, fdm_nhf*, ghostcell*, slice&,bool);
    void force_calc_stl2(lexer*, fdm_nhf*, ghostcell*, slice&,bool);
    void hydrodynamic_viscous_forces_nhflow(lexer*, fdm_nhf*, ghostcell*,slice&, double&,double&,double&,double,double,double,double,double,double,double);
    void force_calc_lsm(lexer*, fdm_nhf*, ghostcell*,slice&);
    void triangulation(lexer*, fdm_nhf*, ghostcell*);
    void reconstruct(lexer*, fdm_nhf*);
    void addpoint(lexer*,fdm_nhf*,int,int);
    void finalize(lexer*,fdm_nhf*);
    double triangle_area(lexer*,double,double,double,double,double,double,double,double,double);
    double clip_edge(double,double);
    double clip_edge_vol(double,double,double);
    bool clip_facet(lexer*,double,double,double,double,double,double,double,double,double, double,double,double,double&,double&,double&,double&);
    void buoyancy_nhflow(lexer*, fdm_nhf*, ghostcell*, double, double&, double&, double&, double&);
    
    nhflow_reinidisc_fsf *pnhfrdisc;
    int *IO,*CR,*CL;
    double *FRK1,*DTT,*LL;
    int *vert, *nflag;
    double *fsf;
};

#endif
