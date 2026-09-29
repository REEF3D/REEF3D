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

#ifndef SFLOW_AMR_SHIP_H_
#define SFLOW_AMR_SHIP_H_

#include"6DOF.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;

using namespace std;

//  Moving body (X 10 2/3) on a refined SFLOW patch: the 6DOF object of the patch.
//  The body itself lives on level 0 (sixdof_sflow); sflow_amr evaluates the body on the
//  patch cells (fs, draft, press) and this class applies it in the patch kernels:
//   X 10 3: surface pressure of the ship, source -WL dp/dx / rho in the momentum equations
//   X 10 2: direct forcing of the body velocity, as sixdof_obj::update_forcing_sflow

class sflow_amr_ship final : public sixdof, public increment
{
public:
    sflow_amr_ship(lexer*);
    virtual ~sflow_amr_ship();

    void start_cfd(lexer*,fdm*,ghostcell*,int,field&,field&,field&,field&,field&,field&,bool) override final {};
    void start_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final {};
    void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool) override final;
    void reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final {};

    void ini(lexer*,ghostcell*) override final {};
    void initialize(lexer*, fdm*, ghostcell*) override final {};
    void initialize(lexer*, fdm2D*, ghostcell*) override final {};
    void initialize(lexer*, fdm_nhf*, ghostcell*) override final {};

    void isource(lexer*,fdm*,ghostcell*) override final {};
    void jsource(lexer*,fdm*,ghostcell*) override final {};
    void ksource(lexer*,fdm*,ghostcell*) override final {};
    void isource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {};
    void jsource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {};
    void ksource(lexer*,fdm_nhf*,ghostcell*,slice&) override final {};

    void isource2D(lexer*,fdm2D*,ghostcell*) override final;
    void jsource2D(lexer*,fdm2D*,ghostcell*) override final;

    double Hsolid(lexer*, slice&);

    slice4 press, draft;

    // body kinematics of level 0: u_fb(0..5), centre of gravity
    double u[6], c[3];
    double zref;

private:
    double alpha[3];
};

#endif
