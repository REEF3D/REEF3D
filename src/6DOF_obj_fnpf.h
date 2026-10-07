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

#ifndef SIXDOF_OBJ_FNPF_H_
#define SIXDOF_OBJ_FNPF_H_

#include"6DOF_obj.h"

//  6DOF body coupled to REEF3D::FNPF: body boundary condition in the sigma grid, hull
//  pressure from Bernoulli with the instantaneous added mass.

class sixdof_obj_fnpf : public sixdof_obj
{
public:
    
    sixdof_obj_fnpf(lexer*, ghostcell*, int);
	virtual ~sixdof_obj_fnpf();
    
    void initialize_fnpf(lexer*, fdm_fnpf*, ghostcell*);
    void solve_eqmotion_fnpf(lexer*, ghostcell*, int, bool);
    void update_position_fnpf(lexer*, ghostcell*, bool);
    void print_fnpf(lexer*, ghostcell*, int);
    void ray_cast_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, slice&);
    void face_data_fnpf(lexer*, fdm_fnpf*, ghostcell*, int, double*, double*, double*, double*);
    void forces_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, double**, bool);
    // forces_fnpf in two parts for several grids (FNPF mesh refinement): the hull triangles
    // whose centroid own(x,y) accepts are integrated on grid (p,c) with the sampling distance
    // del, then the sums are reduced and stored
    struct fnpf_force_sum { double F[3], Mo[3], Am[36], Atot; };
    void forces_fnpf_zero(lexer*, fnpf_force_sum&);
    void forces_fnpf_sum(lexer*, fdm_fnpf*, double*, double**, bool, double, const std::function<bool(double,double)>*, fnpf_force_sum&);
    void forces_fnpf_set(lexer*, ghostcell*, fnpf_force_sum&, bool);
    double fnpf_dsm() const {return DSM;}
    bool fnpf_fixed(lexer*);
    void print_force_fnpf(lexer*);

private:
    
    void externalForces_fnpf(lexer*, ghostcell*, int, bool);
};

#endif
