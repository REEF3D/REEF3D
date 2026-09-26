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

// FNPF interface of the 6DOF kernel: initialisation, RK3 (kernel rk3) or classical RK4
// in step with fnpf_RK3/fnpf_RK4, and the added-mass coupling (M + A) a = F.
// The rigid-body state (p_, c_, h_, e_), objects, mooring and output are shared
// with CFD/NHFLOW; only the geometry (ray_cast_fnpf) and the loads
// (forces_fnpf, 6DOF_obj_fnpf_forces.cpp) are FNPF specific.

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"

void sixdof_obj::initialize_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->mpirank==0)
    {
    cout<<"6DOF_fnpf_ini "<<endl;
    mkdir("./REEF3D_FNPF_6DOF",0777);
    mkdir("./REEF3D_FNPF_6DOF_STL",0777);
    }
    
    // level-set defined objects need the NHFLOW fdm
    if(p->X131>0 || p->X132>0 || p->X133>0 || p->X153>0)
    {
    if(p->mpirank==0)
    cout<<"6DOF FNPF: only triangulated objects (box, cylinder, wedge, hexahedron, STL) are supported"<<endl;
    
    pgc->final();
    exit(1);
    }
    
    if(p->X320>0 && p->mpirank==0)
    cout<<"6DOF FNPF: nets (X 320) are ignored"<<endl;
    
    // mean horizontal cell size
    int num=0;
    DSM=0.0;
    
    SLICELOOP4
    {
    DSM += (p->j_dir==0) ? p->DXN[IP] : 0.5*(p->DXN[IP] + p->DYN[JP]);
    ++num;
    }
    
    DSM = pgc->globalsum(DSM);
    num = pgc->globalisum(num);
    DSM = DSM/double(MAX(num,1));
    
    if(p->X50==1)
	print_ini_vtp(p,pgc);
    
    if(p->X50==2)
    print_ini_stl(p,pgc);
    
    ini_parallel(p,pgc);
	objects_create(p,pgc);
	ini_fbvel(p,pgc);
    
    // mass, CoG and inertia from the surface triangulation
    geometry_parameters_nhflow(p,nullptr,pgc);
    
    iniPosition_RBM(p,pgc);
    quat_matrices(p);
    
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
    
	update_fbvel(p,pgc);
    
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    
	// Mooring (vtk output of the lines is only written for NHFLOW/CFD)
	if(p->X310==0)
	{
		pmooring.push_back(new mooring_void());
	}
	else
	{
		pgc->bcast_int(&p->mooring_count,1);	
		
		Xme.resize(p->mooring_count);
		Yme.resize(p->mooring_count);
		Zme.resize(p->mooring_count);
		Kme.resize(p->mooring_count);
		Mme.resize(p->mooring_count);
		Nme.resize(p->mooring_count);
		
		pmooring.reserve(p->mooring_count);
		
		X311_xen.resize(p->mooring_count,0.0);
		X311_yen.resize(p->mooring_count,0.0);
		X311_zen.resize(p->mooring_count,0.0);

		for (int ii=0; ii<p->mooring_count; ii++)
		{
			if(p->X310==1)
			pmooring.push_back(new mooring_Catenary(ii));
			
			else if(p->X310==2)
			pmooring.push_back(new mooring_barQuasiStatic(ii)); 
			
			else if(p->X310==3)
            pmooring.push_back(new mooring_dynamic(ii));
			
			else if(p->X310==4)
			pmooring.push_back(new mooring_Spring(ii));
			
			Eigen::Vector3d fl(p->X311_xe[ii] - p->xg, p->X311_ye[ii] - p->yg, p->X311_ze[ii] - p->zg);
			fl = R_.transpose()*fl;
			X311_xen[ii] = fl(0);
			X311_yen[ii] = fl(1);
			X311_zen[ii] = fl(2);

			pmooring[ii]->initialize(p,pgc);
		}
	}	
    
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xext=Yext=Zext=Kext=Mext=Next=0.0;
    
    Aadd_.setZero();
    am_on_=false;
    
    for(int s=0; s<3; ++s)
    {
    rk4_p_[s].setZero();
    rk4_c_[s].setZero();
    rk4_h_[s].setZero();
    rk4_e_[s].setZero();
    }
}

bool sixdof_obj::fnpf_fixed(lexer *p)
{
    bool all=true;
    
    for(int n=0; n<6; ++n)
    all = all && p_fixed_dof(p,n);
    
    return all;
}

void sixdof_obj::solve_eqmotion_fnpf(lexer *p, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_fnpf(p,pgc,finalize);
    
    update_forces(p);
    
    // stage-synchronous with fnpf_RK3 (TVD Shu-Osher, as the kernel's rk3) or fnpf_RK4
    if(p->A310==3)
    rk3(p,pgc,iter);
    
    else
    rk4(p,pgc,iter);
}

void sixdof_obj::externalForces_fnpf(lexer *p, ghostcell *pgc, bool finalize)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

	if(p->X310>0)
	mooringForces(p,pgc,1.0);
}

void sixdof_obj::rk4(lexer *p, ghostcell *pgc, int iter)
{
    // classical RK4, stage-synchronous with fnpf_RK4:
    // y_s = y_n + c_s*dt*K_s (c_s = 1/2, 1/2, 1),  y_n+1 = y_n + dt/6*(K1 + 2K2 + 2K3 + K4)
    get_trans(p, pgc, dp_, dc_, p_, c_);    
    get_rot(p, dh_, de_, h_, e_);
    
    if(iter==0)
    {
        pk_ = p_;
        ck_ = c_;
        hk_ = h_;
        ek_ = e_;
    }
    
    if(iter<3)
    {
        rk4_p_[iter] = dp_;
        rk4_c_[iter] = dc_;
        rk4_h_[iter] = dh_;
        rk4_e_[iter] = de_;
        
        const double cs = (iter==2) ? 1.0 : 0.5;
        
        p_ = pk_ + cs*p->dt*dp_;
        c_ = ck_ + cs*p->dt*dc_;
        h_ = hk_ + cs*p->dt*dh_;
        e_ = ek_ + cs*p->dt*de_;
        e_.normalize();
    }
    
    if(iter==3)
    {
        const double w = p->dt/6.0;
        
        p_ = pk_ + w*(rk4_p_[0] + 2.0*rk4_p_[1] + 2.0*rk4_p_[2] + dp_);
        c_ = ck_ + w*(rk4_c_[0] + 2.0*rk4_c_[1] + 2.0*rk4_c_[2] + dc_);
        h_ = hk_ + w*(rk4_h_[0] + 2.0*rk4_h_[1] + 2.0*rk4_h_[2] + dh_);
        e_ = ek_ + w*(rk4_e_[0] + 2.0*rk4_e_[1] + 2.0*rk4_e_[2] + de_);
        e_.normalize();
    }
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
    int flag=1;
    
    if(n==0) flag = p->X11_u;
    if(n==1) flag = p->X11_v;
    if(n==2) flag = p->X11_w;
    if(n==3) flag = p->X11_p;
    if(n==4) flag = p->X11_q;
    if(n==5) flag = p->X11_r;
    
    if(p->j_dir==0 && (n==1 || n==3 || n==5))
    return true;
    
    return (flag!=1);
}

void sixdof_obj::update_position_fnpf(lexer *p, ghostcell *pgc, bool finalize)
{
    quat_matrices(p);
    
    update_Euler_angles(p,pgc);
    
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {
            Eigen::Vector3d point(tri_x0[n][q], tri_y0[n][q], tri_z0[n][q]);
            point = R_*point;
            
            tri_x[n][q] = point(0) + c_(0);
            tri_y[n][q] = point(1) + c_(1);
            tri_z[n][q] = point(2) + c_(2);
        }
	}
    
    omega_B = I_.inverse()*h_;
    omega_I = R_*omega_B;
    
    update_fbvel(p,pgc);
    
    if(p->mpirank==0 && finalize==true)
    {
        cout<<"XG: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
        cout<<"Ue: "<<u_fb(0)<<" Ve: "<< u_fb(1)<<" We: "<< u_fb(2)<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }
}

void sixdof_obj::print_fnpf(lexer *p, ghostcell *pgc, int iter)
{
    saveTimeStep(p,iter);
    
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    
    print_parameter(p,pgc);
}
