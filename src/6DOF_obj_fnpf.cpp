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
#include"ckpt_pack.h"

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
    // index of the initial fluid output (printcount_sixdof); not on a restart from a checkpoint
    // (I 45 1: printed already)
    if(p->I45!=1)
    {
    if(p->X50==1)
    print_vtp(p,pgc);
    
    if(p->X50==2)
    print_stl(p,pgc);
    }
    
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

void sixdof_fluid_fnpf::velocity(int n, const double *xyz, double *uvw)
{
    // each point is interpolated by the rank that owns it (FNPF is decomposed horizontally),
    // the others contribute 0; U, V, W from velcalc_sig of the current stage
    for(int q=0; q<n; ++q)
    {
        const double x = xyz[3*q], y = xyz[3*q+1], z = xyz[3*q+2];
        
        uvw[3*q] = uvw[3*q+1] = uvw[3*q+2] = 0.0;
        
        const bool own = (x>=p->originx && x<p->endx) && (p->j_dir==0 || (y>=p->originy && y<p->endy));
        
        if(own)
        {
            uvw[3*q]   = p->ccipol7V(c->U, c->WL, c->bed, x, y, z);
            uvw[3*q+1] = p->ccipol7V(c->V, c->WL, c->bed, x, y, z);
            uvw[3*q+2] = p->ccipol7V(c->W, c->WL, c->bed, x, y, z);
        }
    }
    
    for(int q=0; q<3*n; ++q)
    uvw[q] = pgc->globalsum(uvw[q]);
}

void sixdof_obj_fnpf::solve_eqmotion_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_fnpf(p,pgc,iter,finalize);
    
    // load models (ship module) sample the FNPF velocity
    sixdof_fluid_fnpf fluid(p,c,pgc);
    pfluid = &fluid;
    
    update_forces(p);
    
    pfluid = nullptr;
    
    // stage-synchronous with fnpf_RK3 (TVD Shu-Osher, as the kernel's rk3) or fnpf_RK4
    if(p->A310==3)
    rk3(p,pgc,iter);
    
    else
    rk4(p,pgc,iter);
}

void sixdof_obj_fnpf::eqmotion_assemble_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, int iter, bool finalize)
{
    externalForces_fnpf(p,pgc,iter,finalize);

    sixdof_fluid_fnpf fluid(p,c,pgc);
    pfluid = &fluid;

    am_coupled_ = true;
    am_assembled_ = false;

    update_forces(p);

    am_coupled_ = false;
    pfluid = nullptr;
}

bool sixdof_obj_fnpf::eqmotion_block_fnpf(Eigen::Matrix<double,6,6> &L, Eigen::Matrix<double,6,1> &r, bool *fixed) const
{
    if(!am_assembled_)
    return false;

    L = Lc_;
    r = rc_;

    for(int n=0; n<6; ++n)
    fixed[n] = fixedc_[n];

    return true;
}

void sixdof_obj_fnpf::eqmotion_advance_fnpf(lexer *p, ghostcell *pgc, int iter, const Eigen::Matrix<double,6,1> *acc)
{
    // a body without a block (no added mass yet, fully fixed) keeps the loads of update_forces
    if(am_assembled_ && acc!=nullptr)
    finish_added_mass(*acc,IIc_,gyroc_,fixedc_);

    am_assembled_ = false;

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
    
    // the pose and velocity: in the 6DOF log files (REEF3D_*_6DOF), no longer on the terminal
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

// checkpoints (P 48, I 45)
void sixdof_obj_fnpf::ckpt_fnpf_state(std::vector<double> &v) const
{
    ckpt_put w{v};
    w.d(2.0);   // version
    w.i(chim_nx); w.i(chim_ny); w.i(chim_count);
    for(double x : {chim_x0,chim_y0,chim_h,chim_bx0,chim_bx1,chim_by0,chim_by1,Fx,Fy,Fz,DSM})
    w.d(x);
}

bool sixdof_obj_fnpf::ckpt_fnpf_restore(const std::vector<double> &v)
{
    ckpt_get r{v};
    double version=0.0;
    r.d(version);
    if(!r.ok || int(version)!=2)
    return false;
    r.i(chim_nx); r.i(chim_ny); r.i(chim_count);
    for(double *x : {&chim_x0,&chim_y0,&chim_h,&chim_bx0,&chim_bx1,&chim_by0,&chim_by1,&Fx,&Fy,&Fz,&DSM})
    r.d(*x);
    return r.ok;
}

void sixdof_obj_fnpf::ckpt_fnpf_rank_state(std::vector<double> &v) const
{
    ckpt_put w{v};
    w.vd(chim_); w.vd(chim2_);
    w.i(int(chim_lev_.size()));
    for(const chim_level &g : chim_lev_)
    {
        w.vd(g.m1); w.vd(g.m2);
        w.i(g.nx); w.i(g.ny);
        w.d(g.x0); w.d(g.y0); w.d(g.h); w.d(g.t);
    }
}

void sixdof_obj_fnpf::ckpt_fnpf_rank_restore(const std::vector<std::vector<double>> &all)
{
    // the means of the ranks (one: as they were); a value per bin: the mean of the ranks that
    // hold one
    std::vector<double> c1, c2;
    std::vector<chim_level> lv;
    std::vector<double> n1;
    std::vector<std::vector<double>> nl;
    
    auto merge = [](std::vector<double> &a, std::vector<double> &n, const std::vector<double> &b)
    {
        if(b.empty())
        return;
        if(a.size()<b.size())
        {
            a.resize(b.size(),0.0);
            n.resize(b.size(),0.0);
        }
        for(size_t q=0; q<b.size(); ++q)
        if(b[q]!=0.0)
        {
            a[q] += b[q];
            n[q] += 1.0;
        }
    };
    
    std::vector<double> nc2;
    std::vector<std::vector<double>> nl2;
    
    for(const std::vector<double> &v : all)
    {
        ckpt_get r{v};
        std::vector<double> a1, a2;
        r.vd(a1); r.vd(a2);
        merge(c1,n1,a1);
        merge(c2,nc2,a2);
        int nlev=0;
        r.i(nlev);
        if(!r.ok || nlev<0)
        continue;
        if((int)lv.size()<nlev)
        {
            lv.resize(nlev);
            nl.resize(nlev);
            nl2.resize(nlev);
        }
        for(int l=0; l<nlev; ++l)
        {
            chim_level g;
            r.vd(g.m1); r.vd(g.m2);
            r.i(g.nx); r.i(g.ny);
            r.d(g.x0); r.d(g.y0); r.d(g.h); r.d(g.t);
            if(!r.ok)
            break;
            if(g.nx>0 && lv[l].nx==0)
            {
                lv[l].nx = g.nx; lv[l].ny = g.ny; lv[l].x0 = g.x0; lv[l].y0 = g.y0; lv[l].h = g.h; lv[l].t = g.t;
            }
            merge(lv[l].m1,nl[l],g.m1);
            merge(lv[l].m2,nl2[l],g.m2);
        }
    }
    
    auto mean = [](std::vector<double> &a, const std::vector<double> &n)
    {
        for(size_t q=0; q<a.size(); ++q)
        if(n[q]>1.5)
        a[q] /= n[q];
    };
    mean(c1,n1);
    mean(c2,nc2);
    for(size_t l=0; l<lv.size(); ++l)
    {
        mean(lv[l].m1,nl[l]);
        mean(lv[l].m2,nl2[l]);
    }
    
    chim_.swap(c1);
    chim2_.swap(c2);
    chim_lev_.swap(lv);
}
