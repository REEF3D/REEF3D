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

#ifndef SIXDOF_OBJ_CFD_H_
#define SIXDOF_OBJ_CFD_H_

#include"6DOF_obj.h"
#include"amr_grid_lookup.h"

class reinidisc_f;

//  Fluid access of the load models (ship module) in REEF3D::CFD: velocity interpolated from the
//  staggered velocity fields
class sixdof_fluid_cfd : public sixdof_fluid
{
public:
    sixdof_fluid_cfd(lexer *pp, fdm *aa, ghostcell *gc) : p(pp), a(aa), pgc(gc) {}
    void velocity(int, const double*, double*) override;
private:
    lexer *p;
    fdm *a;
    ghostcell *pgc;
};

//  6DOF body coupled to REEF3D::CFD: level set of the body on the Cartesian grid, direct
//  forcing and surface force integration (X 60 1 STL triangles, 2 level-set triangulation).

class sixdof_obj_cfd : public sixdof_obj
{
public:
    
    sixdof_obj_cfd(lexer*, ghostcell*, int);
	virtual ~sixdof_obj_cfd();
    
    void solve_eqmotion_cfd(lexer*,fdm*,ghostcell*,int,bool);
    void initialize_cfd(lexer*,fdm*,ghostcell*);
    void update_forcing(lexer*, fdm*, ghostcell*,field&,field&,field&,field&,field&,field&,int);
    // the forcing weights of the solids and the topography at the faces in this stage (a->fbh1-3 of
    // ghostcell::solid_forcing), set by sixdof_cfd and cfd_amr_6dof for update_forcing: the weight of
    // a face of the body is limited so that the weights of the solids and of all bodies at the face
    // add up to 1 at most (a body that reaches into a solid or into another body)
    field *sfh[3] = {nullptr,nullptr,nullptr};
    void hydrodynamic_forces_cfd(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void update_position_3D(lexer*, fdm*, ghostcell*, bool);
    
    // lateral line (X 194): pressure and velocity relative to the hull at the probes, printed with
    // the body output
    void lateral_line(lexer*, fdm*, ghostcell*);
    
    // X 198: fluid-momentum coupling (6DOF_obj_cfd_momentum.cpp): the rigid-body motion fitted to the
    // predicted velocity of the stage in the hull region gives the loads of the body
    void momentum_fit(lexer*, fdm*, ghostcell*, field&, field&, field&, int);
    void momentum_deform(lexer*, fdm*, ghostcell*);
    // ... in parts for the mesh refinement (cfd_amr: the predicted velocity of the patches of the
    // finest level): begin, add the grids, end (sums over the ranks, the solve)
    void momentum_fit_begin(lexer*, int);
    void momentum_fit_add(lexer*, fdm*, field&, field&, field&);
    void momentum_fit_end(lexer*, ghostcell*);
    void momentum_deform_begin();
    void momentum_deform_add(lexer*, fdm*);
    void momentum_deform_end(lexer*, ghostcell*);
    template<class G> void momentum_acc(lexer*, fdm*, G, Eigen::Matrix<double,6,6>&, Eigen::Matrix<double,6,1>&);
    Eigen::Matrix<double,6,1> momentum_solve(lexer*, ghostcell*, Eigen::Matrix<double,6,6>&, Eigen::Matrix<double,6,1>&, bool&);
    Eigen::Matrix<double,6,6> mf_N_ = Eigen::Matrix<double,6,6>::Zero();
    Eigen::Matrix<double,6,1> mf_b_ = Eigen::Matrix<double,6,1>::Zero();
    double mf_c1_ = 0.0, mf_c2_ = 1.0, mf_tn_ = 0.0;
    bool mf_start_ = false, mf_first_ = true;
    void momentum_loads(lexer*, int);
    void momentum_refit(lexer*, fdm*, ghostcell*, field&, field&, field&);
    void momentum_refit_set(lexer*, ghostcell*);
    
    // actuator disks of the load models (body-force propellers): acceleration of the fluid,
    // added to the forcing terms fx, fy, fz of the stage (water part of the cells outside the body)
    void actuator_forcing(lexer*, fdm*, ghostcell*, field&, field&, field&);
    
    using sixdof_obj::print_vtp;

    // mesh refinement (cfd_amr, 6DOF_obj_cfd_amr.cpp): the level set of the body on the grid of a
    // patch, the cells around the patch included, with the workspace of that grid; both work in
    // the box of the body on the grid (box: is, ie, js, je, ks, ke, set by ray_cast_amr)
    void ray_cast_amr(lexer*, fdm*, ghostcell*, fieldint5&, fieldint5&, fieldint5&, int*);
    void reini_amr(lexer*, fdm*, ghostcell*, field&, reinidisc_f*, field4a&, field4a&, field4a&, const int*);
    // ... the loads (X 60 1) from the finest grid at each hull triangle (set by cfd_amr, G 7 1)
    const amr_grid_lookup *amr_loads = nullptr;
    // ... G 7 1: the level-0 cells under the patches (1, per cell of the level-0 arrays): the
    // predicted copy of level 0 does not force the body there (set by cfd_amr for its step)
    const char *amr_cover = nullptr;
    // ... G 7 1: stage iter of the body with the step of a level (p->dt): the loads of its last
    // stage, the pose and the body velocities (the hull is cast on the patches, not on level 0)
    void amr_stage_cfd(lexer*, fdm*, ghostcell*, int);

private:
    
    void ini_parameter_stl(lexer*, fdm*, ghostcell*);
    void externalForces_cfd(lexer*, fdm*, ghostcell*, double, bool);
    void netForces_cfd(lexer*, fdm*, ghostcell*, double, bool);
    double Hsolidface(lexer*, fdm*, int,int,int);
    void face_normal(lexer*, fdm*, int, double&, double&, double&);
    double Hsolidface_t(lexer*, fdm*, int,int,int);
    void geometry_parameters(lexer*, fdm*, ghostcell*);
    void geometry_ls(lexer*, fdm*, ghostcell*);
    void forces_stl(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void forces_lsm(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void triangulation(lexer*, fdm*, ghostcell*, field&);
    void reconstruct(lexer*, fdm*, field&);
    void addpoint(lexer*,fdm*,int,int);
    void forces_lsm_calc(lexer* p, fdm *a, ghostcell *pgc,int,bool);
    void print_force(lexer*,fdm*,ghostcell*);
    void print_ini(lexer*,fdm*,ghostcell*);
    void update_trimesh_3D(lexer*, fdm*, ghostcell*, bool);
    void ray_cast(lexer*, fdm*, ghostcell*);
    void reini_RK2(lexer*, fdm*, ghostcell*, field&);
    
    fieldint5 cutl,cutr,fbio;
    reinidisc *prdisc;
    field4a f, frk1, L, dt;
    fieldint5 vertice, nodeflag;
    field5 eta;
    
    void actuator_point(lexer*, int, Eigen::Vector3d&, double&);
    double actuator_fluid_fraction(lexer*, fdm*, int);
    bool actuator_warned = false;
};

#endif
