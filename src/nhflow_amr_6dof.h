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

#ifndef NHFLOW_AMR_6DOF_H_
#define NHFLOW_AMR_6DOF_H_

#include"6DOF.h"
#include"increment.h"

class sixdof_nhflow;

using namespace std;

//  The floating bodies (6DOF_nhflow, X 10 1/2) on a refined NHFLOW grid (nhflow_amr): the
//  sixdof object of a patch.  The body is advanced on level 0 (sixdof_nhflow, before the patches
//  take their forcing; with subcycling, G 7 1, by the finest level: nhflow_amr_sub.cpp); a patch
//  casts the hull on its own sigma grid (level set FB, ray-cast workspace of the patch size) and adds the direct forcing of the rigid-body velocity, as
//  sixdof_nhflow::start_twoway/reforce_nhflow do it on level 0.  The loads are integrated once,
//  on level 0, with every hull triangle sampled on the finest grid (sixdof_obj::amr_grid_nhflow).
//  Without bodies (nullptr) all calls do nothing.

class nhflow_amr_6dof final : public sixdof, public increment
{
public:
    nhflow_amr_6dof(lexer*, sixdof_nhflow*, double);
    virtual ~nhflow_amr_6dof();

    // the hull on the grid: level set FB, solid flags DF and DFF as nhflow_forcing::forcing
    void body(lexer*, fdm_nhf*, ghostcell*);

    void start_cfd(lexer*,fdm*,ghostcell*,int,field&,field&,field&,field&,field&,field&,bool) override final {}
    void start_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final;
    void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool) override final {}

    void reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final;

    void ini(lexer*,ghostcell*) override final {}
    void initialize(lexer*, fdm*, ghostcell*) override final {}
    void initialize(lexer*, fdm2D*, ghostcell*) override final {}
    void initialize(lexer*, fdm_nhf*, ghostcell*) override final {}

    void isource(lexer*,fdm*,ghostcell*) override final {}
    void jsource(lexer*,fdm*,ghostcell*) override final {}
    void ksource(lexer*,fdm*,ghostcell*) override final {}

    void isource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {}
    void jsource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {}
    void ksource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {}

    void isource2D(lexer*,fdm2D*,ghostcell*) override final {}
    void jsource2D(lexer*,fdm2D*,ghostcell*) override final {}

private:
    void ray_cast(lexer*, fdm_nhf*, ghostcell*, int);

    sixdof_nhflow *b;
    lexer *pp;
    int *IO, *CL, *CR;      // ray-cast workspace of the grid
    double dsm;             // mean cell size of the grid (band of the level set)
    int n7;
};

#endif
