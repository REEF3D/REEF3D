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
    void hydrodynamic_forces_cfd(lexer*, fdm*, ghostcell*,field&,field&,field&,int,bool);
    void update_position_3D(lexer*, fdm*, ghostcell*, bool);
    
    // actuator disks of the load models (body-force propellers): acceleration of the fluid,
    // added to the forcing terms fx, fy, fz of the stage (water part of the cells outside the body)
    void actuator_forcing(lexer*, fdm*, ghostcell*, field&, field&, field&);
    
    using sixdof_obj::print_vtp;

private:
    
    void ini_parameter_stl(lexer*, fdm*, ghostcell*);
    void externalForces_cfd(lexer*, fdm*, ghostcell*, double, bool);
    void netForces_cfd(lexer*, fdm*, ghostcell*, double, bool);
    double Hsolidface(lexer*, fdm*, int,int,int);
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
    
    fieldint4 cutl,cutr,fbio;
    reinidisc *prdisc;
    field4a f, frk1, L, dt;
    fieldint4 vertice, nodeflag;
    field4 eta;
    
    void actuator_point(lexer*, int, Eigen::Vector3d&, double&);
    double actuator_fluid_fraction(lexer*, fdm*, int);
    bool actuator_warned = false;
};

#endif
