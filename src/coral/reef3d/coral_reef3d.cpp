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

#include"coral/reef3d/coral_reef3d.h"
#include"coral/reef3d/coral_fluid_sampler.h"
#include"coral/io/coral_coupling.h"
#include"coral/io/coral_input.h"
#include"coral/core/coral_version.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<stdexcept>
#include<sys/stat.h>

namespace
{
    std::ostream nullout(nullptr);      // messages of the ranks other than 0
}

coral_reef3d::coral_reef3d(lexer *p, ghostcell *pgc) : p_(p), pgc_(pgc)
{
    // rank 0 reads coral.txt, all ranks build the same model
    std::string text;
    int ok = 1;
    if(p->mpirank==0)
    {
        std::ifstream in("coral.txt");
        if(!in)
            ok = 0;
        else
        {
            std::stringstream ss;
            ss << in.rdbuf();
            text = ss.str();
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, pgc->mpi_comm);
    if(!ok)
        fail("Z 40: coral.txt not found in the case directory");
    int n = (int)text.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, pgc->mpi_comm);
    text.resize(n);
    if(n > 0)
        MPI_Bcast(&text[0], n, MPI_CHAR, 0, pgc->mpi_comm);

    case_ = std::make_unique<coral::coral_case>();
    try
    {
        case_->read_string(text, ".");
    }
    catch(const std::exception &e)
    {
        fail(e.what());
    }
    std::ostream &log = p->mpirank==0 ? std::cout : nullout;
    cp_ = std::make_unique<coral::coupling>(*case_, log, p->mpirank==0);
    log<<coral::banner()<<": "<<case_->m.nbody()<<" bodies, "<<case_->m.ncon()<<" constraint equations, "
       <<case_->m.load_list.size()<<" loads, max. step "<<case_->run.dt<<" s"<<std::endl;
    // fluid bodies and members: CORAL is the external dynamics of the 6DOF bodies of the fluid solver
    if(case_->fluid || !case_->fluid_members.empty())
    {
        sixdof_external::attach(this);
        if(case_->fluid)
        {
            log<<"CORAL: "<<case_->fluid_bodies.size()<<" fluid bod"<<(case_->fluid_bodies.size()>1 ? "ies" : "y")<<" (6DOF bodies";
            for(const auto &fb : case_->fluid_bodies)
                log<<" "<<fb.n<<" = "<<case_->m.body(fb.body).name<<(fb.one_way ? " one-way" : "");
            log<<"), loads of order "<<case_->fluid_order<<(case_->fluid_predictor ? ", predictor" : "")<<std::endl;
        }
        if(case_->run.static_ini)
            log<<"CORAL: static equilibrium with fluid bodies: the fluid loads are not known yet (zero)"<<std::endl;
        const coral::Vec3 g(p->W20, p->W21, p->W22);
        if((case_->m.gravity - g).norm() > 1.0e-6*std::max(1.0, g.norm()))
            log<<"CORAL: WARNING gravity "<<case_->m.gravity.transpose()<<" in coral.txt, "<<g.transpose()
               <<" in ctrl.txt (W 20-22); the fluid loads come without the weight"<<std::endl;
    }
}

coral_reef3d::~coral_reef3d()
{
    if(sixdof_external::instance()==this)
        sixdof_external::attach(nullptr);
}

void coral_reef3d::fail(const std::string &msg)
{
    std::cerr<<"\n!!! REEF3D::CORAL: "<<msg<<" !!!\n"<<std::endl;
    MPI_Abort(pgc_->mpi_comm, 1);
}

void coral_reef3d::set_sampler(coral_fluid_sampler *s)
{
    sampler_.reset(s);
}

void coral_reef3d::sample_kinematics(double t)
{
    if(!cp_->has_kinematics())
        return;
    if(!sampler_)
    {
        if(!sampler_told_ && p_->mpirank==0)
            std::cout<<"CORAL: WARNING fluid_kinematics needs FNPF or NHFLOW, the waves of coral.txt are used"<<std::endl;
        sampler_told_ = true;
        return;
    }
    std::vector<coral::Vec3> pts;
    cp_->kinematics_points(pts);
    const int n = (int)pts.size();
    if(n == 0)
        return;
    std::vector<double> xyz(3*n), loc, sum(5*n, 0.0);
    for(int q=0; q<n; ++q)
        for(int d=0; d<3; ++d)
            xyz[3*q+d] = pts[q](d);
    sampler_->sample(xyz, loc);
    MPI_Allreduce(loc.data(), sum.data(), 5*n, MPI_DOUBLE, MPI_SUM, pgc_->mpi_comm);
    bcast(sum.data(), 5*n);
    std::vector<coral::Vec3> u(n, coral::Vec3::Zero());
    std::vector<double> eta(n, 0.0);
    std::vector<int> in(n, 0);
    int nin = 0;
    for(int q=0; q<n; ++q)
    {
        const double c = sum[5*q+4];
        if(c > 0.5)
        {
            u[q] = coral::Vec3(sum[5*q], sum[5*q+1], sum[5*q+2])/c;
            eta[q] = sum[5*q+3]/c;
            in[q] = 1;
            ++nin;
        }
    }
    cp_->kinematics(t, u, eta, in);
    if(!sampler_told_ && p_->mpirank==0)
    {
        std::cout<<"CORAL: fluid kinematics from "<<sampler_->name()<<" at "<<n<<" points ("<<nin<<" in the fluid domain)"<<std::endl;
        sampler_told_ = true;
    }
}

void coral_reef3d::ini(lexer *p, ghostcell*)
{
    if(p->mpirank==0)
        mkdir("./REEF3D_CORAL", 0777);
    try
    {
        cp_->record("./REEF3D_CORAL");
        sample_kinematics(p->simtime);
        cp_->ini(p->simtime, "./REEF3D_CORAL");
    }
    catch(const std::exception &e)
    {
        fail(e.what());
    }
}

void coral_reef3d::step(lexer *p, ghostcell*)
{
    const double t1 = p->simtime + p->dt;
    // the fluid bodies have already advanced CORAL within the fluid step
    if(t1 <= cp_->time() + 1.0e-12*std::max(1.0, std::fabs(t1)))
        return;
    try
    {
        // the fluid at the end of its step: the CORAL step interpolates the kinematics
        sample_kinematics(t1);
        cp_->step_to(t1);
        check_ranks("step");
    }
    catch(const std::exception &e)
    {
        fail(e.what());
    }
}

// ---- external dynamics of the 6DOF bodies (6DOF_external.h)

bool coral_reef3d::drives(int n) const
{
    return cp_->drives(n);
}

bool coral_reef3d::one_way(int n) const
{
    return cp_->one_way(n);
}

bool coral_reef3d::predictor() const
{
    return cp_->predictor();
}

bool coral_reef3d::strips(int n, Eigen::Vector3d &a, Eigen::Vector3d &e, double &L, int &ns) const
{
    return cp_->strips(n, a, e, L, ns);
}

void coral_reef3d::fluid_strips(int n, double t, const std::vector<double> &F)
{
    const int m = (int)F.size();
    std::vector<double> sum(m, 0.0);
    MPI_Allreduce(F.data(), sum.data(), m, MPI_DOUBLE, MPI_SUM, pgc_->mpi_comm);
    bcast(sum.data(), m);
    coral::VecX Fv(m);
    for(int i=0; i<m; ++i)
        Fv(i) = sum[i];
    cp_->fluid_strips(n, t, Fv);
}

void coral_reef3d::connect(int n, double mass, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Matrix3d &I)
{
    // mass properties of rank 0 on all ranks
    double v[22];
    v[0] = mass;
    for(int i=0; i<3; ++i)
        v[1+i] = c(i);
    for(int i=0; i<9; ++i)
    {
        v[4+i] = R(i/3, i%3);
        v[13+i] = I(i/3, i%3);
    }
    bcast(v, 22);
    Eigen::Vector3d cb(v[1], v[2], v[3]);
    Eigen::Matrix3d Rb, Ib;
    for(int i=0; i<9; ++i)
    {
        Rb(i/3, i%3) = v[4+i];
        Ib(i/3, i%3) = v[13+i];
    }
    cp_->connect(n, v[0], cb, Rb, Ib);
}

void coral_reef3d::fluid_loads(double t, const std::vector<vec6> &F, const std::vector<mat6> &A, bool provisional)
{
    // the loads (summed over the ranks by the fluid solver) and the added mass of rank 0 on all ranks
    std::vector<vec6> Fb(F);
    std::vector<mat6> Ab(A);
    if(p_->mpi_size > 1)
    {
        std::vector<double> v;
        v.reserve(6*Fb.size() + 36*Ab.size() + 1);
        v.push_back(t);
        for(const auto &f : Fb)
            for(int i=0; i<6; ++i)
                v.push_back(f(i));
        for(const auto &a : Ab)
            for(int i=0; i<36; ++i)
                v.push_back(a(i/6, i%6));
        bcast(v.data(), (int)v.size());
        t = v[0];
        int k = 1;
        for(auto &f : Fb)
            for(int i=0; i<6; ++i)
                f(i) = v[k++];
        for(auto &a : Ab)
            for(int i=0; i<36; ++i)
                a(i/6, i%6) = v[k++];
    }
    try
    {
        // the kinematics of the same fluid state first
        if(!provisional)
            sample_kinematics(t);
        cp_->fluid_loads(t, Fb, Ab, provisional);
    }
    catch(const std::exception &e)
    {
        fail(e.what());
    }
}

void coral_reef3d::advance(double t0, double t1)
{
    try
    {
        cp_->advance(t0, t1);
        check_ranks("advance");
    }
    catch(const std::exception &e)
    {
        fail(e.what());
    }
}

void coral_reef3d::predict(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w)
{
    cp_->predict(n, t, c, R, v, w);
}

void coral_reef3d::state(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w)
{
    cp_->state(n, t, c, R, v, w);
}

void coral_reef3d::bcast(double *v, int n)
{
    if(p_->mpi_size > 1 && n > 0)
        MPI_Bcast(v, n, MPI_DOUBLE, 0, pgc_->mpi_comm);
}

void coral_reef3d::check_ranks(const char *where)
{
    if(p_->mpi_size < 2)
        return;
    // pose and velocities of the bodies and the generalized coordinates: the same on every rank
    const coral::model &m = case_->m;
    std::vector<double> v;
    v.push_back(cp_->time());
    for(const auto &b : m.bodies)
    {
        for(int i=0; i<3; ++i)
            v.push_back(b.x(i));
        for(int i=0; i<9; ++i)
            v.push_back(b.R(i/3, i%3));
        for(int i=0; i<3; ++i)
            v.push_back(b.v(i));
        for(int i=0; i<3; ++i)
            v.push_back(b.w(i));
    }
    for(int i=0; i<m.cq.size(); ++i)
        v.push_back(m.cq(i));
    for(int i=0; i<m.cv.size(); ++i)
        v.push_back(m.cv(i));
    const int n = (int)v.size();
    std::vector<double> lo(n), hi(n);
    MPI_Allreduce(v.data(), lo.data(), n, MPI_DOUBLE, MPI_MIN, pgc_->mpi_comm);
    MPI_Allreduce(v.data(), hi.data(), n, MPI_DOUBLE, MPI_MAX, pgc_->mpi_comm);
    double d = 0.0;
    for(int i=0; i<n; ++i)
        d = std::max(d, hi[i] - lo[i]);
    if(d > 0.0 && !rank_diff_told_ && p_->mpirank==0)
        std::cout<<"CORAL: WARNING the ranks differ in the state of the model after "<<where<<" at t = "<<cp_->time()
                 <<" (largest difference "<<d<<")"<<std::endl;
    if(d > 0.0)
        rank_diff_told_ = true;
}
