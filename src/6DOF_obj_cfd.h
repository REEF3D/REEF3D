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
    void print_vtp(lexer*,fdm*,ghostcell*);
    void pvtp(lexer*,int);
    void update_trimesh_3D(lexer*, fdm*, ghostcell*, bool);
    void ray_cast(lexer*, fdm*, ghostcell*);
    void reini_RK2(lexer*, fdm*, ghostcell*, field&);
    
    fieldint5 cutl,cutr,fbio;
    reinidisc *prdisc;
    field4a f, frk1, L, dt;
    fieldint5 vertice, nodeflag;
    field5 eta;
};

#endif
