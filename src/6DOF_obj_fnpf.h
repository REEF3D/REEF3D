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
#include<vector>

//  Fluid access of the load models (ship module) in REEF3D::FNPF: velocity of the potential flow
//  (grad phi on the sigma grid)
class sixdof_fluid_fnpf : public sixdof_fluid
{
public:
    sixdof_fluid_fnpf(lexer *pp, fdm_fnpf *cc, ghostcell *gc) : p(pp), c(cc), pgc(gc) {}
    void velocity(int, const double*, double*) override;
private:
    lexer *p;
    fdm_fnpf *c;
    ghostcell *pgc;
};

//  6DOF body coupled to REEF3D::FNPF: body boundary condition in the sigma grid, hull
//  pressure from Bernoulli with the instantaneous added mass.

class sixdof_obj_fnpf : public sixdof_obj
{
public:
    
    sixdof_obj_fnpf(lexer*, ghostcell*, int);
	virtual ~sixdof_obj_fnpf();
    
    void initialize_fnpf(lexer*, fdm_fnpf*, ghostcell*);
    struct fnpf_force_sum { double F[3], Mo[3], Am[36], Atot; };
    void solve_eqmotion_fnpf(lexer*, fdm_fnpf*, ghostcell*, int, bool);
    // several bodies (X 20 > 1): the stage in two parts around the coupled added-mass solve
    // of fnpf_6DOF::motion: loads and the body's block, then the accelerations and the RK stage
    void eqmotion_assemble_fnpf(lexer*, fdm_fnpf*, ghostcell*, int, bool);
    bool eqmotion_block_fnpf(Eigen::Matrix<double,6,6>&, Eigen::Matrix<double,6,1>&, bool*) const;
    void eqmotion_advance_fnpf(lexer*, ghostcell*, int, const Eigen::Matrix<double,6,1>*);
    // cross added mass: -rho int psi_m N dS over this hull for the unit modes psi_m of another body
    void added_mass_block_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, double**, Eigen::Matrix<double,6,6>&);
    // the same from an integral already summed on this rank (several grids): reduce and store
    void added_mass_from_sum(ghostcell*, fnpf_force_sum&, Eigen::Matrix<double,6,6>&);
    bool dof_fixed_fnpf(lexer *p, int n) {return p_fixed_dof(p,n);}
    void update_position_fnpf(lexer*, ghostcell*, bool);
    // mesh refinement with subcycling (fnpf_amr, G 7 1): the stage of the predicted copy, the
    // loads frozen (6DOF_obj_amr.cpp)
    void amr_predict_fnpf(lexer*, ghostcell*, int);
    void print_fnpf(lexer*, ghostcell*, int);
    void ray_cast_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, slice&);
    void face_data_fnpf(lexer*, fdm_fnpf*, ghostcell*, int, double*, double*, double*, double*);
    void forces_fnpf(lexer*, fdm_fnpf*, ghostcell*, double*, double**, bool);
    // body-following time derivative chi = phi_t + V.grad(phi) (V = u_c + w x r): free-surface
    // Dirichlet data of the psi_0 solve, chi = phi_t + V.grad(phi) outside the footprint
    void chi_fsf(lexer*, fdm_fnpf*, slice&, slice&);
    // X 18 tau: remove the running mean (time constant tau) of chi at body-fixed points from the
    // free-surface data next to the body (the time mean of the exact chi, the derivative of phi
    // following a body point, vanishes in any statistically steady state)
    void chi_mean(lexer*, ghostcell*, slice&, slice&);
    // X 18 on the patches of FNPF mesh refinement (one body-frame grid per refinement level, cell
    // size of that level; the patch lexer, the level)
    void chi_mean_patch(lexer*, ghostcell*, slice&, slice&, int);
    static bool chi_on;   // psi_0 is chi (one body), otherwise phi_t with the m-terms neglected
    // forces_fnpf in two parts for several grids (FNPF mesh refinement): the hull triangles
    // whose centroid own(x,y) accepts are integrated on grid (p,c) with the sampling distance
    // del, then the sums are reduced and stored
    void forces_fnpf_zero(lexer*, fnpf_force_sum&);
    void forces_fnpf_sum(lexer*, fdm_fnpf*, double*, double**, bool, double, const std::function<bool(double,double)>*, fnpf_force_sum&);
    void forces_fnpf_set(lexer*, ghostcell*, fnpf_force_sum&, bool);
    double fnpf_dsm() const {return DSM;}
    bool fnpf_fixed(lexer*);
    void print_force_fnpf(lexer*);
    // checkpoints (P 48, I 45): the FNPF part of the state (the grid of the running mean of chi of
    // X 18, the hull loads of the last psi solve, the mean cell size), after that of the kernel
    // (sixdof_obj::ckpt_state)
    void ckpt_fnpf_state(std::vector<double>&) const;
    bool ckpt_fnpf_restore(const std::vector<double>&);
    // the running means of X 18 are kept by every rank for its columns (rank blob); on another
    // decomposition the mean over the ranks that hold a value
    void ckpt_fnpf_rank_state(std::vector<double>&) const;
    void ckpt_fnpf_rank_restore(const std::vector<std::vector<double>>&);

private:
    
    void externalForces_fnpf(lexer*, ghostcell*, int, bool);
    
    // X 18: running mean of chi on a body-frame grid of cell size chim_h around the body
    std::vector<double> chim_, chim2_;
    int chim_nx=0, chim_ny=0, chim_count=-1;
    double chim_x0=0.0, chim_y0=0.0, chim_h=1.0;
    double chim_bx0=0.0, chim_bx1=0.0, chim_by0=0.0, chim_by1=0.0;   // body extent in the body frame
    struct chim_level
    {
        std::vector<double> m1, m2;
        int nx=0, ny=0;
        double x0=0.0, y0=0.0, h=1.0, t=-1.0e20;
    };
    std::vector<chim_level> chim_lev_;   // X 18 on the patches, by refinement level
};

#endif
