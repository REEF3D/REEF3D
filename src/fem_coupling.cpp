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

#include"fem_coupling.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#include<mpi.h>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdio>
#include<algorithm>
#include<sys/stat.h>
#include<sys/types.h>

fem_coupling::fem_coupling(lexer *p, ghostcell *pgc) : surf_version(-1), dxmin(0.0), rho_w(1000.0), initialised(false), force_scale(1.0), nstep(0), printtime(0.0), printcount(0), starttime(0.0)
{
    // rank 0 reads fem.dat, everybody parses the broadcast content
    std::string content;
    int len = 0;

    if(p->mpirank==0)
    {
        std::ifstream f("fem.dat");
        if(!f)
        len = -1;
        else
        {
            std::stringstream ss;
            ss<<f.rdbuf();
            content = ss.str();
            len = (int)content.size();
        }
    }

    MPI_Bcast(&len,1,MPI_INT,0,pgc->mpi_comm);

    if(len<0)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! Z 30: fem.dat not found in the case directory !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    content.resize(len);
    if(len>0)
    MPI_Bcast(&content[0],len,MPI_CHAR,0,pgc->mpi_comm);

    rho_w = p->W1;

    // smallest fluid cell size: spacing of the Lagrangian points, default element size
    double h = 1.0e20;
    for(int ii=0; ii<p->knox; ++ii)
    h = std::min(h,p->DXN[ii+marge]);
    if(p->j_dir==1)
    for(int jj=0; jj<p->knoy; ++jj)
    h = std::min(h,p->DYN[jj+marge]);
    for(int kk=0; kk<p->knoz; ++kk)
    h = std::min(h,p->DZN[kk+marge]);
    dxmin = pgc->globalmin(h);

    try
    {
        fs.set_gravity(fem_solid::Vec3(p->W20,p->W21,p->W22));
        if(p->j_dir==0)
        fs.set_plane_strain(true);

        std::istringstream is(content);
        fs.read(is);

        // element size from the fluid grid ('resolution'), one layer across a 2D slice
        if(!fs.lattice_given())
        fs.set_default_spacing(fs.coupling().resolution*dxmin, p->j_dir==0 ? p->global_ymax-p->global_ymin : -1.0);

        fs.build();
    }
    catch(std::exception& e)
    {
        if(p->mpirank==0)
        std::cout<<"\n!!! "<<e.what()<<" !!!\n"<<std::endl;
        MPI_Abort(pgc->mpi_comm,1);
    }

    force_scale = (p->j_dir==0) ? 1.0/fs.lattice_h(1) : 1.0;
    outdir = "./REEF3D_FEM";

    if(p->mpirank==0)
    {
        mkdir(outdir.c_str(),0777);
        fs.info(std::cout);
    }
}

void fem_coupling::first_call(lexer *p, fdm *a, ghostcell *pgc)
{
    initialised = true;

    // supports on the bed and the solids of the fluid grid
    // the boundaries of the fluid domain are walls for free bodies and debris
    if(fs.coupling().walls)
    {
        fs.add_contact_plane(fem_solid::Vec3(1,0,0), p->global_xmin);
        fs.add_contact_plane(fem_solid::Vec3(-1,0,0),-p->global_xmax);
        fs.add_contact_plane(fem_solid::Vec3(0,0,1), p->global_zmin);
        if(p->j_dir==1)
        {
            fs.add_contact_plane(fem_solid::Vec3(0,1,0), p->global_ymin);
            fs.add_contact_plane(fem_solid::Vec3(0,-1,0),-p->global_ymax);
        }
    }

    if(fs.coupling().fix_bed)
    {
        std::vector<double> lv(fs.nnode(),1.0e30);
        for(int i=0; i<fs.nnode(); ++i)
        {
            fem_solid::Vec3 x = fs.ref_pos(i);
            if(p->j_dir==0) x(1) = p->YP[marge];
            if(x(0)>=p->originx && x(0)<p->endx && (p->j_dir==0 || (x(1)>=p->originy && x(1)<p->endy)) && x(2)>=p->originz && x(2)<p->endz)
            lv[i] = std::min(p->ccipol4a(a->topo,x(0),x(1),x(2)),p->ccipol4a(a->solid,x(0),x(1),x(2)));
        }
        if(!lv.empty())
        MPI_Allreduce(MPI_IN_PLACE,lv.data(),(int)lv.size(),MPI_DOUBLE,MPI_MIN,pgc->mpi_comm);

        std::vector<int> bed;
        const double tol = 0.3*fs.hmin();
        for(int i=0; i<fs.nnode(); ++i)
        if(lv[i]<=tol && !fs.node_rigid_material(i))
        bed.push_back(i);
        fs.fix_nodes(bed);
        if(p->mpirank==0)
        std::cout<<"FEM: fix bed: "<<bed.size()<<" nodes on the bed / solids of the grid"<<std::endl;
    }

    std::vector<std::string> warn;
    warnings(p,warn);

    // check mode: report and stop
    if(fs.coupling().check)
    {
        fem_solid::check_info ci = fs.check();
        for(const std::string& w : warn) ci.warnings.push_back(w);

        // wetted surface at the start
        double wet = 0.0;
        for(const fem_solid::face& f : fs.surface())
        {
            fem_solid::Vec3 c = 0.25*(fs.pos(f.n[0])+fs.pos(f.n[1])+fs.pos(f.n[2])+fs.pos(f.n[3]));
            if(p->j_dir==0) c(1) = p->YP[marge];
            if(c(0)>=p->originx && c(0)<p->endx && (p->j_dir==0 || (c(1)>=p->originy && c(1)<p->endy)) && c(2)>=p->originz && c(2)<p->endz)
            if(p->ccipol4(a->phi,c(0),c(1),c(2))>=0.0)
            wet += 0.5*((fs.pos(f.n[2])-fs.pos(f.n[0])).cross(fs.pos(f.n[3])-fs.pos(f.n[1]))).norm();
        }
        wet = pgc->globalsum(wet);

        if(p->mpirank==0)
        {
            std::ostringstream os;
            fs.write_check(os,ci);
            os<<"  wetted surface at the start: "<<wet<<" m2\n";
            os<<"  fluid: smallest cell "<<dxmin<<" m, element size / fluid cell "<<fs.hmin()/dxmin<<"\n";
            if(p->dt>0.0)
            os<<"  substeps per fluid step of "<<p->dt<<" s: about "<<std::max(1,(int)std::ceil(p->dt/(0.5*fs.critical_dt())))<<"\n";
            std::cout<<"\n"<<os.str()<<std::endl;
            std::ofstream f((outdir+"/REEF3D_FEM_check.txt").c_str());
            f<<os.str();
            fs.write_vtu(outdir+"/REEF3D-FEM-check.vtu");
            std::cout<<"FEM: check written to "<<outdir<<"/REEF3D_FEM_check.txt and REEF3D-FEM-check.vtu, stopping (remove 'check' from fem.dat to run)"<<std::endl;
        }
        pgc->final(false);
    }

    if(p->mpirank==0)
    for(const std::string& w : warn)
    std::cout<<"FEM WARNING: "<<w<<std::endl;

    ini_points(p,pgc);

    // settling before the flow: self weight and the pressure of the initial
    // fluid (a structure standing in still water starts in equilibrium)
    if(fs.coupling().settle && (fs.has_supports() || fs.ground()))
    {
        std::vector<fem_solid::Vec3> F0;
        // initial pressure field, or hydrostatic below the initial free
        // surface if the flow solver starts without one (I 12 0)
        pressure_loads(p,a,pgc,F0,p->I12<1);
        fem_solid::Vec3 Ft = fem_solid::Vec3::Zero();
        fs.clear_loads();
        for(int i=0; i<fs.nnode(); ++i)
        {
            fs.add_load(i,F0[i]);
            Ft += F0[i];
        }

        double res = 0.0;
        const bool ok = fs.settle(200000,1.0e-4,&res,!fs.ground(),true);
        const fem_solid::Vec3 R = fs.support_force();
        fs.update_utilisation();
        const double u0 = fs.max_utilisation();
        fs.clear_loads();
        if(p->mpirank==0 && u0>1.0)
        std::cout<<"FEM WARNING: the structure is overloaded already in the initial state (utilisation "<<u0
                 <<" under self weight and the still water): it will crack or yield at the start of the run"<<std::endl;

        if(p->mpirank==0)
        {
            double mass = 0.0;
            for(int i=0; i<fs.nnode(); ++i) mass += fs.mass(i);
            std::cout<<"FEM: settled under self weight"<<(Ft.norm()>1.0e-3*mass*9.81 ? " and the initial fluid pressure" : "")
                     <<(ok ? "" : " (NOT converged)")<<": max displacement "<<fs.max_displacement()*1000.0
                     <<" mm, weight "<<mass*9.81/1000.0<<" kN, initial fluid force "<<Ft(0)/1000.0<<" "<<Ft(1)/1000.0<<" "<<Ft(2)/1000.0
                     <<" kN, support force "<<R(0)/1000.0<<" "<<R(1)/1000.0<<" "<<R(2)/1000.0<<" kN"<<std::endl;
        }
    }

    // structural damping at the first natural frequency
    fs.prepare_damping();
    if(p->mpirank==0 && fs.damping_alpha()>0.0)
    std::cout<<"FEM: structural damping "<<100.0*fs.damping_ratio()<<" % at the first natural frequency "<<fs.damping_frequency()<<" Hz (dry)"<<std::endl;

    if(p->mpirank==0)
    {
        std::cout<<"FEM: "<<pts.size()<<" Lagrangian points on the surface, critical solid time step "<<fs.critical_dt()<<" s"<<std::endl;

        std::ofstream lg((outdir+"/REEF3D_FEM_log.dat").c_str());
        lg<<"# time  substeps  intact_elements  eroded  debris  F_fluid_x F_fluid_y F_fluid_z  F_support_x F_support_y F_support_z  max_vonMises  max_displacement  kinetic_energy  dissipated_energy  (SI units, support force = force of the structure on its supports"<<(p->j_dir==0 ? ", 2D: forces of the FEM slice" : "")<<")\n";

        for(const fem_solid::monitor& mo : fs.monitors())
        {
            std::ofstream ts((outdir+"/REEF3D_FEM_monitor_"+mo.name+".dat").c_str());
            const fem_solid::Vec3& X = fs.ref_pos(mo.node);
            ts<<"# node at "<<X(0)<<" "<<X(1)<<" "<<X(2)<<"\n# time  dx dy dz  vx vy vz\n";
        }
    }

    print(p);
}

void fem_coupling::warnings(lexer *p, std::vector<std::string>& w)
{
    std::ostringstream os;

    double rmin = 1.0e30;
    for(int k=0; k<fs.material_count(); ++k)
    if(!fs.mat(k).rigid)
    rmin = std::min(rmin,fs.mat(k).rho);
    if(rmin < 1.1*rho_w)
    {
        os.str(""); os<<"a deformable material is lighter than 1.1 x water ("<<rmin<<" kg/m3): light deformable structures can become unstable; "
                       "for floating debris use a rigid material ('material timber C24 rigid' or 'material rigid 450')";
        w.push_back(os.str());
    }
    if(fs.rigid_material_deformable())
    w.push_back("a rigid material is used in a body with supports or together with deformable materials: it is computed as elastic there");
    if(fs.hmin() > 2.01*dxmin)
    {
        os.str(""); os<<"elements ("<<fs.hmin()<<" m) are more than twice the fluid cells ("<<dxmin<<" m): the flow is resolved, but stresses are coarse";
        w.push_back(os.str());
    }
    if(fs.hmin() < 0.24*dxmin)
    {
        os.str(""); os<<"elements ("<<fs.hmin()<<" m) are much finer than the fluid cells ("<<dxmin<<" m): many substeps, details smaller than a fluid cell feel no flow";
        w.push_back(os.str());
    }

    fem_solid::Vec3 lo = fem_solid::Vec3::Constant(1.0e300), hi = fem_solid::Vec3::Constant(-1.0e300);
    for(int i=0; i<fs.nnode(); ++i)
    {
        lo = lo.cwiseMin(fs.ref_pos(i));
        hi = hi.cwiseMax(fs.ref_pos(i));
    }
    const double tol = 1.0e-6;
    if(lo(0)<p->global_xmin-tol || hi(0)>p->global_xmax+tol || lo(2)<p->global_zmin-tol || hi(2)>p->global_zmax+tol
       || (p->j_dir==1 && (lo(1)<p->global_ymin-tol || hi(1)>p->global_ymax+tol)))
    w.push_back("the structure reaches outside the fluid domain: parts outside get no fluid loads");

    if(!fs.has_supports() && !fs.ground() && fs.n_rigid()==0)
    w.push_back("the structure has no supports ('fix base', 'fix bed' or a fix box) and no ground: it will fall or drift");
}

void fem_coupling::update_summary(lexer *p, bool write)
{
    const double t = fs.time();
    const fem_solid::Vec3 R = fs.support_force();
    const fem_solid::Vec3 M = fs.support_moment();
    const fem_solid::Vec3 F = fs.total_load();

    const double shear = std::sqrt(R(0)*R(0)+R(1)*R(1))*force_scale;
    const double moment = std::sqrt(M(0)*M(0)+M(1)*M(1))*force_scale;
    const double fluid = F.norm()*force_scale;
    const double disp = fs.max_displacement();

    if(shear>sm.shear) {sm.shear = shear; sm.t_shear = t;}
    if(moment>sm.moment) {sm.moment = moment; sm.t_moment = t;}
    if(fluid>sm.fluid) {sm.fluid = fluid; sm.t_fluid = t;}
    if(disp>sm.disp) {sm.disp = disp; sm.t_disp = t;}

    bool event = false;

    if(write || nstep%5==0)
    {
        fs.update_utilisation();
        int e = -1;
        const double u = fs.max_utilisation(&e);
        if(e>=0 && u>sm.util) {sm.util = u; sm.t_util = t; sm.x_util = fs.elem_centre(e);}

        if(sm.t_crack<0.0)
        {
            int ed = -1, ep = -1;
            const double d = fs.max_damage(&ed);
            const double pl = fs.max_plastic_strain(&ep);
            if(d>0.05 || pl>1.0e-4)
            {
                sm.t_crack = t;
                sm.yielding = !(d>0.05);
                sm.x_crack = fs.elem_centre(d>0.05 ? ed : ep);
                event = true;
                if(p->mpirank==0)
                std::cout<<"FEM: first "<<(sm.yielding ? "yielding" : "cracking")<<" at t = "<<t<<" s at ("<<sm.x_crack.transpose()<<")"<<std::endl;
            }
        }
    }

    if(sm.t_fail<0.0 && fs.n_eroded()>0)
    {
        sm.t_fail = t;
        for(int e=0; e<fs.nelem(); ++e)
        if(!fs.elem(e).alive) {sm.x_fail = fs.elem_centre(e); break;}
        event = true;
        if(p->mpirank==0)
        std::cout<<"FEM: first element failure at t = "<<t<<" s at ("<<sm.x_fail.transpose()<<")"<<std::endl;
    }

    if((write || event || nstep%20==0) && p->mpirank==0)
    write_summary(p);
}

void fem_coupling::write_summary(lexer *p)
{
    const double ef = fs.eroded_mass_fraction();
    const char* unit = (p->j_dir==0) ? " per metre width" : "";

    std::ostringstream st;
    if(ef>0.5) st<<"COLLAPSED: "<<std::setprecision(3)<<100.0*ef<<" % of the mass has failed";
    else if(ef>0.0) st<<"PARTLY FAILED: "<<std::setprecision(3)<<100.0*ef<<" % of the mass has failed";
    else if(sm.t_crack>=0.0) st<<(sm.yielding ? "YIELDED, no failure" : "CRACKED, no failure");
    else if(sm.util>=0.0) st<<"INTACT, max utilisation "<<std::setprecision(3)<<sm.util;
    else st<<"INTACT (elastic material, max von Mises "<<std::setprecision(4)<<fs.max_vonmises()/1.0e6<<" MPa)";

    std::ofstream f((outdir+"/REEF3D_FEM_summary.txt").c_str());
    f<<std::setprecision(4);
    f<<"REEF3D FEM summary at t = "<<fs.time()<<" s\n\n";

    if(fs.n_deformable()>0 || fs.n_rigid()==0)
    {
    f<<"status:                  "<<st.str()<<"\n";
    f<<"max base shear:          "<<sm.shear/1000.0<<" kN"<<unit<<"  at t = "<<sm.t_shear<<" s\n";
    f<<"max overturning moment:  "<<sm.moment/1000.0<<" kNm"<<unit<<"  at t = "<<sm.t_moment<<" s  (about the base centre "<<fs.base_centre().transpose()<<")\n";
    f<<"max total fluid force:   "<<sm.fluid/1000.0<<" kN"<<unit<<"  at t = "<<sm.t_fluid<<" s\n";
    f<<"max displacement:        "<<sm.disp*1000.0<<" mm  at t = "<<sm.t_disp<<" s\n";
    if(sm.util>=0.0)
    f<<"max utilisation:         "<<sm.util<<"  at t = "<<sm.t_util<<" s at ("<<sm.x_util.transpose()<<")   (stress / strength, > 1: cracking or yielding)\n";
    else
    f<<"max utilisation:         - (elastic materials only)\n";
    if(sm.t_crack>=0.0)
    f<<"first "<<(sm.yielding ? "yielding:          " : "cracking:          ")<<"t = "<<sm.t_crack<<" s at ("<<sm.x_crack.transpose()<<")\n";
    else
    f<<"first cracking:          none\n";
    if(sm.t_fail>=0.0)
    f<<"first element failure:   t = "<<sm.t_fail<<" s at ("<<sm.x_fail.transpose()<<")\n";
    else
    f<<"first element failure:   none\n";
    f<<"failed mass:             "<<100.0*ef<<" %,  debris particles "<<fs.debris().size()<<"\n";
    }

    // rigid bodies: floating and drifting debris
    if(fs.n_rigid()>0)
    {
        f<<(fs.n_deformable()>0 ? "\n" : "")<<"rigid bodies (debris): "<<fs.n_rigid()<<(p->j_dir==0 ? "   (2D: mass and forces per metre width)" : "")<<"\n";
        f<<"  body   mass [kg]  density   start centre                 current centre               max speed [m/s]  max drift [m]  max contact force [kN]\n";
        const int nmax = std::min(fs.n_rigid(),50);
        for(int k=0; k<nmax; ++k)
        {
            const fem_solid::rigid_body& rb = fs.rigid(k);
            f<<"  "<<std::setw(4)<<k+1<<"  "<<std::setw(10)<<rb.M*force_scale<<"  "<<std::setw(7)<<rb.M/std::max(rb.Vol,1.0e-30)
             <<"   ("<<std::setw(7)<<rb.c0(0)<<" "<<std::setw(7)<<rb.c0(1)<<" "<<std::setw(7)<<rb.c0(2)<<")"
             <<"   ("<<std::setw(7)<<rb.c(0)<<" "<<std::setw(7)<<rb.c(1)<<" "<<std::setw(7)<<rb.c(2)<<")"
             <<"   "<<std::setw(10)<<rb.vmax<<"   "<<std::setw(10)<<rb.dmax<<"   "<<std::setw(10)<<rb.fcmax*force_scale/1000.0<<"\n";
        }
        if(fs.n_rigid()>nmax)
        f<<"  ... "<<fs.n_rigid()-nmax<<" more\n";
    }
}

fem_coupling::~fem_coupling()
{
}

void fem_coupling::ini_points(lexer *p, ghostcell *pgc)
{
    // Lagrangian points: every surface face is split into n1 x n2 points with
    // a spacing below the smallest fluid cell. In 2D (j_dir 0) the faces
    // normal to y are the out-of-plane boundaries and carry no points.
    (void)pgc;
    pts.clear();

    const std::vector<fem_solid::face>& F = fs.surface();
    const int npf = fs.coupling().points_per_face;

    for(int f=0; f<(int)F.size(); ++f)
    {
        const fem_solid::Vec3& x0 = fs.ref_pos(F[f].n[0]);
        const fem_solid::Vec3& x1 = fs.ref_pos(F[f].n[1]);
        const fem_solid::Vec3& x3 = fs.ref_pos(F[f].n[3]);

        const fem_solid::Vec3 nrm = (fs.ref_pos(F[f].n[2])-x0).cross(x3-x1);
        if(p->j_dir==0 && std::fabs(nrm(1))>0.5*nrm.norm())
        continue;

        int n1 = npf>0 ? npf : std::max(1,(int)std::ceil((x1-x0).norm()/(0.9*dxmin)));
        int n2 = npf>0 ? npf : std::max(1,(int)std::ceil((x3-x0).norm()/(0.9*dxmin)));

        // 2D: one point across the out-of-plane direction
        if(p->j_dir==0)
        {
            const fem_solid::Vec3 e1 = x1-x0, e2 = x3-x0;
            if(std::fabs(e1(1))>0.5*e1.norm()) n1 = 1;
            if(std::fabs(e2(1))>0.5*e2.norm()) n2 = 1;
        }

        for(int a=0; a<n1; ++a)
        for(int b=0; b<n2; ++b)
        {
            lpoint q;
            q.face = f;
            q.s = (a+0.5)/double(n1);
            q.t = (b+0.5)/double(n2);
            q.frac = 1.0/double(n1*n2);
            pts.push_back(q);
        }
    }

    surf_version = fs.surface_version();

    buf.assign(14*(pts.size()+fs.debris().size()) + 2*size_t(fs.nnode()),0.0);
    fdeb.assign(fs.debris().size(),fem_solid::Vec3::Zero());
}

void fem_coupling::point_state(int q, fem_solid::Vec3& xp, fem_solid::Vec3& vp, fem_solid::Vec3& n, double& A) const
{
    const lpoint& L = pts[q];
    const fem_solid::face& F = fs.surface()[L.face];

    // node order n0 -> n1 -> n2 -> n3 around the face, s along n0-n1, t along n0-n3
    const double w[4] = {(1.0-L.s)*(1.0-L.t), L.s*(1.0-L.t), L.s*L.t, (1.0-L.s)*L.t};

    xp.setZero();
    vp.setZero();
    for(int a=0; a<4; ++a)
    {
        xp += w[a]*fs.pos(F.n[a]);
        vp += w[a]*fs.vel_mean(F.n[a]);
    }

    const fem_solid::Vec3 c = (fs.pos(F.n[2])-fs.pos(F.n[0])).cross(fs.pos(F.n[3])-fs.pos(F.n[1]));
    const double cn = c.norm();
    n = cn>0.0 ? fem_solid::Vec3(c/cn) : fem_solid::Vec3::Zero();
    A = 0.5*cn*L.frac;
}

double fem_coupling::kernel(double r) const
{
    // 3-point kernel of Roma et al. (1999), as in the FSI strips
    r = std::fabs(r);
    if(r<=0.5)
    return (1.0 + std::sqrt(1.0 - 3.0*r*r))/3.0;
    if(r<=1.5)
    return (5.0 - 3.0*r - std::sqrt(1.0 - 3.0*(1.0-r)*(1.0-r)))/6.0;
    return 0.0;
}

void fem_coupling::print(lexer *p)
{
    if(p->mpirank!=0)
    return;

    const double dtp = fs.coupling().print_dt>0.0 ? fs.coupling().print_dt : p->Z31;

    if(dtp>0.0 && fs.time()>=printtime-1.0e-12)
    {
        char name[256];
        std::snprintf(name,sizeof(name),"%s/REEF3D-FEM-%08d.vtu",outdir.c_str(),printcount);
        try
        {
            fs.write_vtu(name);
        }
        catch(std::exception& e)
        {
            std::cout<<e.what()<<std::endl;
        }
        ++printcount;
        printtime += dtp;
        write_summary(p);
    }
}
