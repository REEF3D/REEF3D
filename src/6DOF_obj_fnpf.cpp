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
#include"6DOF_obj_fnpf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include"mooring_void.h"
#include"mooring_barQuasiStatic.h"
#include"mooring_Catenary.h"
#include"mooring_Spring.h"
#include"mooring_dynamic.h"

sixdof_obj_fnpf::sixdof_obj_fnpf(lexer *p, ghostcell *pgc, int number) : sixdof_obj(p,pgc,number)
{
}

sixdof_obj_fnpf::~sixdof_obj_fnpf()
{
}


void sixdof_obj_fnpf::initialize_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->mpirank==0)
    {
    cout<<"6DOF_fnpf_ini "<<endl;
    mkdir("./REEF3D_FNPF_6DOF",0777);
    mkdir("./REEF3D_FNPF_6DOF_STL",0777);
    
    if(p->X310>0)
    mkdir("./REEF3D_FNPF_6DOF_Mooring",0777);
    }
    
    // all objects are triangulated; mass, CoG and inertia from the surface triangles
    
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
    geometry_parameters_surface(p,pgc,0.0);
    
    iniPosition_RBM(p,pgc);
    quat_matrices(p);
    
    rb.update_omega();
    
	update_fbvel(p,pgc);
    
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
    
    // t = 0 frame of the body after the mooring set-up, so both series share the
    // index of the initial fluid output (printcount_sixdof)
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xext=Yext=Zext=Kext=Mext=Next=0.0;
    
    // power take-off (X 500)
    ini_pto(p,pgc);

    Aadd_.setZero();
    am_on_=false;
    
}

bool sixdof_obj_fnpf::fnpf_fixed(lexer *p)
{
    bool all=true;
    
    for(int n=0; n<6; ++n)
    all = all && p_fixed_dof(p,n);
    
    return all;
}

void sixdof_obj_fnpf::solve_eqmotion_fnpf(lexer *p, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_fnpf(p,pgc,iter,finalize);
    
    update_forces(p);
    
    // stage-synchronous with fnpf_RK3 (TVD Shu-Osher, as the kernel's rk3) or fnpf_RK4
    if(p->A310==3)
    rk3(p,pgc,iter);
    
    else
    rk4(p,pgc,iter);
}

void sixdof_obj_fnpf::externalForces_fnpf(lexer *p, ghostcell *pgc, int iter, bool finalize)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

	if(p->X310>0)
	mooringForces(p,pgc,1.0);

    if(pto_on_)
    pto_forces(p,pgc,iter);
}

void sixdof_obj_fnpf::update_position_fnpf(lexer *p, ghostcell *pgc, bool finalize)
{
    quat_matrices(p);
    
    rb.euler_angles();
    
    geom.transform(R_,c_);
    
    rb.update_omega();
    
    update_fbvel(p,pgc);
    
    if(p->mpirank==0 && finalize==true)
    {
        cout<<"XG: "<<c_(0)<<" YG: "<<c_(1)<<" ZG: "<<c_(2)<<" phi: "<<phi*(180.0/PI)<<" theta: "<<theta*(180.0/PI)<<" psi: "<<psi*(180.0/PI)<<endl;
        cout<<"Ue: "<<u_fb(0)<<" Ve: "<< u_fb(1)<<" We: "<< u_fb(2)<<" Pe: "<<omega_I(0)<<" Qe: "<<omega_I(1)<<" Re: "<<omega_I(2)<<endl;
    }
}

void sixdof_obj_fnpf::print_fnpf(lexer *p, ghostcell *pgc, int iter)
{
    saveTimeStep(p,iter);
    
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    
    print_parameter(p,pgc);
}

void sixdof_obj_fnpf::print_force_fnpf(lexer *p)
{
    // same columns as CFD/NHFLOW; potential flow: pressure part = total, no viscous part
    if(p->mpirank==0)
    printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
    <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Fx<<" \t "<<Fy<<" \t "<<Fz<<" \t "<<0.0<<" \t "<<0.0<<" \t "<<0.0<<endl;
}
