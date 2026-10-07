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

#ifndef FNPF_BODY_H_
#define FNPF_BODY_H_

class lexer;
class fdm_fnpf;
class ghostcell;
class solver;
class fnpf_laplace;
class fnpf_fsf;
class slice;
class fnpf_amr;
class sixdof_obj;

#include<vector>

// Interface between the FNPF time stepping and resolved bodies.
// The base class is the "no body" case: every hook is a no-op and the Laplace solver
// is returned unchanged, so fnpf_RK3/fnpf_RK4 call the hooks unconditionally and
// contain no body logic. create() returns fnpf_6DOF for X 10 > 0.
//
// initialize() is called once at t = 0 (inidisc_step2, after the sigma grid is set up),
// so the body output starts together with the initial fluid output.
// Per RK stage the scheme calls
//   stage()   with the tendencies deta/dt, dFifsf/dt of the current state, before the
//             stage value is formed (body loads and body RK stage),
//   surface() after the stage values eta/Fifsf are formed, before fsfdisc/sigma_update,
//   velocity() after the velocities of the end of the step are computed.
// Geometry and body-band extrapolation around the Laplace solve are hidden behind the
// solver returned by laplace().

class fnpf_body
{
public:
    static fnpf_body* create(lexer*, fdm_fnpf*, ghostcell*);
    
    virtual ~fnpf_body(){}
    
    // once at t = 0 after the initial sigma grid, before the initial output
    virtual void initialize(lexer*, fdm_fnpf*, ghostcell*){}
    
    virtual void stage(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_fsf*, slice&, slice&, int){}
    virtual void surface(lexer*, fdm_fnpf*, ghostcell*, slice&, slice&, int, int){}
    virtual fnpf_laplace* laplace(fnpf_laplace *plap){return plap;}
    // after velcalc_sig of the end of a time step (output, time step): velocities at the
    // body nodes
    virtual void velocity(lexer*, fdm_fnpf*, ghostcell*){}
    
    // mesh refinement (fnpf_amr, G 1): the body is represented on every grid of the
    // hierarchy; fnpf_amr calls these around its composite phi solve and after each patch
    // stage value; the psi solves of the loads run on all grids as well
    virtual bool present() const {return false;}
    virtual void amr_attach(fnpf_amr*){}
    virtual void amr_grids(lexer*, ghostcell*){}
    virtual void amr_geometry(lexer*, fdm_fnpf*, ghostcell*){}
    virtual void amr_post_solve(lexer*, fdm_fnpf*, ghostcell*, double*){}
    virtual void amr_surface(lexer*, ghostcell*, int, slice&, slice&){}
    virtual void amr_bodies(std::vector<sixdof_obj*>&){}
};

#endif
