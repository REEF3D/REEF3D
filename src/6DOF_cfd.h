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

#ifndef SIXDOF_CFD_H_
#define SIXDOF_CFD_H_

#include"6DOF.h"
#include"6DOF_obj_cfd.h"
#include<vector>

class mooring;
class ddweno_f_nug;

using namespace std;

class sixdof_cfd final : public sixdof, public increment
{
public:
    // runtime coupling (reefcouple, I 243): state of the bodies
    int body_count() const override {return int(fb_obj.size());}
    void body_state(int, double*) const override;
    
	sixdof_cfd(lexer*, fdm*, ghostcell*);
	virtual ~sixdof_cfd();

    void start_cfd(lexer*,fdm*,ghostcell*,int,field&,field&,field&,field&,field&,field&,bool) override final;
    void start_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final;
    void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool) override final;
    
    void reforce_nhflow(lexer*,fdm_nhf*,ghostcell*,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool) override final {};
    
    void ckpt_state(lexer*, std::vector<double>&) override final;
    bool ckpt_restore(lexer*, ghostcell*, const std::vector<double>&) override final;
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

    // the bodies, for the CFD mesh refinement (cfd_amr)
    int objects() const {return number6DOF;}
    sixdof_obj_cfd* object(int nb) {return fb_obj[nb];}
    // ... G 7 1: level 0 steps with a predicted copy of the bodies (the loads of the last stage of
    // the finest level, no output); the finest level advances them (sixdof_obj_cfd::amr_stage_cfd)
    bool amr_predict = false;

private:
   void setup(lexer*,fdm*,ghostcell*);
   
    // external solver of the bodies (6DOF_external.h, REEF3D::CORAL fluid_body; 6DOF_cfd_external.cpp)
    void ext_connect(lexer*, ghostcell*);
    void ext_nodes(lexer*, int, double&, double&, int&) const;
    void ext_loads(lexer*, double, bool);
    bool ext_drives(int) const;
    void ext_stage(lexer*, int);
    class sixdof_external *pext = nullptr;
    bool ext_checked = false;
    bool ext_twoway = false;
   
    int number6DOF;
    vector<sixdof_obj_cfd*> fb_obj;
    
    // several bodies: the level set of all bodies (minimum) in a->fb after their update; the ray cast
    // of a body overwrites a->fb with its own, the forcing of the body uses that one
    void fb_combine(lexer*, fdm*, int);
    void fb_finish(lexer*, fdm*, ghostcell*);
    field4a *fball = nullptr;

};

#endif
