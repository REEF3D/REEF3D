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

#ifndef FNPF_6DOF_H_
#define FNPF_6DOF_H_

#include"increment.h"
#include"slice4.h"
#include<vector>

class lexer;
class fdm_fnpf;
class ghostcell;
class solver;
class fnpf_laplace;
class fnpf_fsf;
class fnpf_bed_update;
class sixdof_obj;
class slice;

using namespace std;

// Resolved rigid bodies in REEF3D::FNPF (X 10 1), driven by the shared 6DOF kernel.
//
// Geometry: vertical ray cast of the trimesh per sigma column -> body nodes (FBF),
// staircase Neumann faces in the Laplace matrix, body footprint where the
// free-surface node lies inside the body. In the footprint eta and Fifsf are a
// smooth (harmonic) extension of the surrounding free surface, so the sigma grid
// stays regular; they carry no physics there.
//
// Loads: psi = phi_t is solved with the same operator as phi,
//   psi_0: Dirichlet phi_t at the free surface (from the RK tendencies), Neumann
//          w x (w x r).n on the hull (m-terms neglected)
//   psi_j: Dirichlet 0, Neumann N_j (unit modes) -> added mass A_ij
// and the body is advanced with (M + A) a = F_0 + F_ext, stage-synchronous with
// fnpf_RK3 or fnpf_RK4. A is refreshed once per time step (first stage).

class fnpf_6DOF : public increment
{
public:
    fnpf_6DOF(lexer*, fdm_fnpf*, ghostcell*);
    virtual ~fnpf_6DOF();
    
    void ini(lexer*, fdm_fnpf*, ghostcell*);
    
    // at the start of RK stage iter, from the current state and its tendencies
    void forces(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_laplace*, fnpf_fsf*, slice&, slice&, int);
    void motion(lexer*, fdm_fnpf*, ghostcell*, int);
    
    // after the stage values eta/Fifsf are formed
    void footprint(lexer*, fdm_fnpf*, ghostcell*, slice&, slice&, int, int);
    
    // after sigma_update, before the Laplace solve
    void geometry(lexer*, fdm_fnpf*, ghostcell*);
    
    // after a Laplace solve: fill the body band for the cross terms and sampling
    void extrapolate(lexer*, fdm_fnpf*, ghostcell*, double*);
    
    bool initialized;
    
private:
    void solve_psi(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_laplace*, fnpf_fsf*, double*, slice&);
    void zero_face(lexer*, fdm_fnpf*);
    void exchange_face(lexer*, fdm_fnpf*, ghostcell*);
    
    vector<sixdof_obj*> fb_obj;
    int nbody;
    int gcval;
    int footcount;
    
    double *psi0;
    double *psi[6];
    double *zero;
    int *mark;
    
    slice4 foot,psiD,zeroslice;
    slice4 eta_ext,fi_ext;
    bool ext_ini;
    
    fnpf_bed_update *pbed;
};

#endif
