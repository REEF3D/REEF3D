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
Architect: Hans Bihs
--------------------------------------------------------------------*/


#include"6DOF_obj.h"
#include"6DOF_obj_fnpf.h"
#include"6DOF_deform.h"
#include"lexer.h"
#include"ghostcell.h"

//  Mesh refinement with subcycling (G 7 1: nhflow_amr, fnpf_amr): the finest level advances the
//  body; the coarser levels step before it with a predicted copy of it.  amr_save keeps the
//  rigid-body state, amr_restore puts it back with the trimesh and the body velocities.

void sixdof_obj::amr_save()
{
    rb_amr = rb;
    ext_on_amr_ = ext_on_;
    ext_v_amr_ = ext_v_;
    ext_w_amr_ = ext_w_;
    am_prev_amr_ = am_prev_;
    
    if(fish_!=nullptr)
    fish_t_amr_ = fish_->time();
    
    if(dc_on_)
    {
    const sixdof_deform::correction &c = fish_->corr();
    const double v[6] = {c.tx, c.tl, c.psi, c.vx, c.vl, c.psit};
    for(int q=0; q<6; ++q)
    {
    dc_amr_[q] = v[q];
    dc_amr_[6+q] = dcn_[q];
    }
    for(int q=0; q<3; ++q)
    {
    dc_amr_[12+q] = dck_[q];
    dc_amr_[15+q] = dcr_[q];
    }
    }
}

void sixdof_obj::pose_update(lexer *p, ghostcell *pgc)
{
    rb.euler_angles();
    geom.transform(R_,c_);
    rb.update_omega();
    update_fbvel(p,pgc);
}

void sixdof_obj::amr_restore(lexer *p, ghostcell *pgc)
{
    rb = rb_amr;
    ext_on_ = ext_on_amr_;
    ext_v_ = ext_v_amr_;
    ext_w_ = ext_w_amr_;
    am_prev_ = am_prev_amr_;
    
    if(fish_!=nullptr)
    deform_time(fish_t_amr_);
    
    if(dc_on_)
    {
    sixdof_deform::correction c;
    c.tx = dc_amr_[0]; c.tl = dc_amr_[1]; c.psi = dc_amr_[2];
    c.vx = dc_amr_[3]; c.vl = dc_amr_[4]; c.psit = dc_amr_[5];
    fish_->set_corr(c);
    for(int q=0; q<6; ++q)
    dcn_[q] = dc_amr_[6+q];
    for(int q=0; q<3; ++q)
    {
    dck_[q] = dc_amr_[12+q];
    dcr_[q] = dc_amr_[15+q];
    }
    }
    
    geom.transform(R_,c_);
    update_fbvel(p,pgc);
}

// FNPF: the RK stage of the predicted copy with the loads and the added mass of the last solve
void sixdof_obj_fnpf::amr_predict_fnpf(lexer *p, ghostcell *pgc, int iter)
{
    if(p->A310==3)
    rk3(p,pgc,iter);
    
    else
    rk4(p,pgc,iter);
    
    update_position_fnpf(p,pgc,false);
}
