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

#ifndef SIXDOF_NHFLOW_H_
#define SIXDOF_NHFLOW_H_

#include"6DOF.h"
#include<vector>
#include"increment.h"
#include"6DOF_obj_nhflow.h"

class lexer;
class fdm2D;
class fdm_nhf;
class ghostcell;
class slice;
class fdm;

using namespace std;

class sixdof_nhflow final : public sixdof, public increment
{
public:
	
    sixdof_nhflow(lexer*, ghostcell*);
    virtual ~sixdof_nhflow();
    
    void start_cfd(lexer*,fdm*,ghostcell*,int,field&,field&,field&,field&,field&,field&,bool) override final;
    void start_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final;
    void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool) override final;
    
    void reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final;
    
    // impermeable membranes (X 330) attached to the floating body
    void membrane_forcing_nhflow(lexer*,fdm_nhf*,ghostcell*,double,double*,double*,double*,slice&) override final;
    void membrane_reaction_nhflow(lexer*,fdm_nhf*,ghostcell*,double,slice&,bool) override final;
    bool membrane_iterated() override final;
    void membrane_reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,double,double*,double*,double*,slice&) override final;
    bool membrane_couple_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double,slice&,int) override final;
    
    void start_twoway(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,slice&,slice&,bool);
    void start_oneway(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,slice&,slice&,bool);
    void start_wavemaker(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,slice&,slice&,bool);
    void start_shipwave(lexer*,fdm_nhf*,ghostcell*,int,bool);
       
	void ini(lexer*,ghostcell*) override final;
    void initialize(lexer*, fdm*, ghostcell*) override final;
    void initialize(lexer*, fdm2D*, ghostcell*) override final;
    void initialize(lexer*, fdm_nhf*, ghostcell*) override final;
	
    
    void isource(lexer*,fdm*,ghostcell*) override final;
    void jsource(lexer*,fdm*,ghostcell*) override final;
    void ksource(lexer*,fdm*,ghostcell*) override final;
    
    void isource(lexer*,fdm_nhf*,ghostcell*,slice&) override final;
    void jsource(lexer*,fdm_nhf*,ghostcell*,slice&) override final;
    void ksource(lexer*,fdm_nhf*,ghostcell*,slice&) override final;
    
    void isource2D(lexer*,fdm2D*,ghostcell*) override final;
    void jsource2D(lexer*,fdm2D*,ghostcell*) override final;
    
    // the bodies, for the NHFLOW mesh refinement (nhflow_amr)
    int objects() const {return number6DOF;}
    sixdof_obj_nhflow* object(int nb) {return fb_obj[nb];}
    
    // subcycling (G 7 1): level 0 steps with a predicted copy of the bodies (X 10 1/2): the RK
    // stages with the loads frozen, the trimesh, ray cast and forcing; no loads, no output
    bool amr_predict = false;
    
private:
	
    // hires gradient
    double limiter(double v1, double v2);
    double starttime;
    
    double denom,val,r,phival;
    double dfdx_plus,dfdx_min,dfdy_plus,dfdy_min,dfdx,dfdy;


    int number6DOF;
    vector<sixdof_obj_nhflow*> fb_obj;

    slice4 press;


};

#endif
