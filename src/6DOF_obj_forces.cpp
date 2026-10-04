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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"

void sixdof_obj_cfd::hydrodynamic_forces_cfd(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter, bool finalize)
{
    if(p->X60==1)
    forces_stl(p,a,pgc,uvel,vvel,wvel,iter,finalize);
    
    if(p->X60==2)
    forces_lsm(p,a,pgc,uvel,vvel,wvel,iter,finalize);
}

void sixdof_obj::update_forces(lexer *p)
{
    // Forces in inertial system: external loads, linear damping, DOF modes
    double Fext[6] = {Xext + Xe, Yext + Ye, Zext + Ze, Kext + Ke, Mext + Me, Next + Ne};
    
    // load models (ship module): evaluated with the state of the stage
    if(!pload.empty())
    {
        double Fl[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        
        for(size_t ql=0; ql<pload.size(); ++ql)
        pload[ql]->add_load(p,rb,geom,Fl);
        
        for(int qn=0; qn<6; ++qn)
        Fext[qn] += Fl[qn];
    }
    
    rb.assemble_loads(Fext);
    
    if(Ffb_(0)!=Ffb_(0))
    cout<<"Ffb_(0)....###"<<endl;
    
    if(Ffb_(1)!=Ffb_(1))
    cout<<"Ffb_(1)....###"<<endl;
    
    if(Ffb_(2)!=Ffb_(2))
    cout<<"Ffb_(2)....###"<<endl;
    
    
    if(Mfb_(0)!=Mfb_(0))
    cout<<"Mfb_(0)....###"<<endl;
    
    if(Mfb_(1)!=Mfb_(1))
    cout<<"Mfb_(1)....###"<<endl;
    
    if(Mfb_(2)!=Mfb_(2))
    cout<<"Mfb_(2)....###"<<endl;
    
    // FNPF: instantaneous added mass (and implicit PTO terms, X 500 2) on the left-hand side
    if(am_on_ || (pto_on_ && pto_implicit_))
    apply_added_mass(p);
}

void sixdof_obj::apply_added_mass(lexer *p)
{
    // Rigid body with the instantaneous added mass A (inertial frame, moments about the CoG):
    //   [M I + A_tt   A_tr    ] [a    ]   [ F                ]
    //   [A_rt      I_I + A_rr ] [alpha] = [ Mo - w x (I_I w) ]
    // update_forces has already assembled F and Mo (hydrodynamic without the
    // acceleration part, gravity, mooring, damping). The kernel integrates
    // dp/dt = F and dh_B/dt = 2 Gdot G^T h + R^T Mo, so handing it
    //   F* = M a,   Mo* = I_I alpha + w x (I_I w)
    // gives dh_B/dt = I_B alpha_B and leaves get_trans/get_rot untouched.
    
    const Eigen::Matrix3d II = R_*I_*R_.transpose();
    const Eigen::Vector3d w = omega_I;
    const Eigen::Vector3d gyro = w.cross(II*w);
    
    Eigen::Matrix<double,6,6> L = Aadd_;
    L.block<3,3>(0,0) += Mass_fb*Eigen::Matrix3d::Identity();
    L.block<3,3>(3,3) += II;
    
    Eigen::Matrix<double,6,1> r;
    r.head<3>() = Ffb_;
    r.tail<3>() = Mfb_ - gyro;
    
    // PTO Jacobians (X 500 2)
    pto_implicit(p,L,r);

    bool fixed[6];
    
    for(int n=0; n<6; ++n)
    {
    fixed[n] = p_fixed_dof(p,n);
    
    if(fixed[n])
    {
    L.row(n).setZero();
    L.col(n).setZero();
    L(n,n) = 1.0;
    r(n) = 0.0;
    }
    }
    
    const Eigen::Matrix<double,6,1> acc = L.partialPivLu().solve(r);
    
    Ffb_ = Mass_fb*acc.head<3>();
    Mfb_ = II*acc.tail<3>() + gyro;
    
    // fixed DOFs keep the kernel's convention of zero load
    for(int n=0; n<3; ++n)
    {
    if(fixed[n])
    Ffb_(n) = 0.0;
    
    if(fixed[n+3])
    Mfb_(n) = 0.0;
    }
}

bool sixdof_obj::p_fixed_dof(lexer *p, int n)
{
    // free DOF: X11 flag 1 (2 = prescribed via motionext, 0 = fixed)
    // 2D: sway, roll and yaw do not exist
    return rb.fixed(n);
}

