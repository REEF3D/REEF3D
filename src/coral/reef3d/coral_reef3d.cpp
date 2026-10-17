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
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_static.h"
#include"coral/linear/coral_modes.h"
#include"coral/io/coral_modes_out.h"
#include"coral/core/coral_version.h"
#include"coral/io/coral_input.h"
#include"coral/io/coral_output.h"
#include"coral/hydro/coral_fluid_load.h"
#include"coral/core/coral_rotation.h"
#include"lexer.h"
#include"ghostcell.h"
#include<mpi.h>
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<sys/stat.h>

coral_reef3d::coral_reef3d(lexer *p, ghostcell *pgc)
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
        fail(pgc, "Z 40: coral.txt not found in the case directory");
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
        fail(pgc, e.what());
    }
    ti_ = std::make_unique<coral::timeint>(case_->run.ti);
    p_ = p;
    pgc_ = pgc;
    if(p->mpirank==0)
    {
        std::cout<<coral::banner()<<": "<<case_->m.nbody()<<" bodies, "<<case_->m.ncon()<<" constraint equations, "
                 <<case_->m.load_list.size()<<" loads, max. step "<<case_->run.dt<<" s"<<std::endl;
    }
    // fluid bodies: CORAL is the external dynamics of the 6DOF bodies of the fluid solver
    if(case_->fluid)
    {
        std::vector<int> b;
        for(const auto &fb : case_->fluid_bodies)
            b.push_back(fb.body);
        track_.set_bodies(b);
        sixdof_external::attach(this);
        if(p->mpirank==0)
        {
            std::cout<<"CORAL: "<<b.size()<<" fluid bod"<<(b.size()>1 ? "ies" : "y")<<" (6DOF bodies";
            for(const auto &fb : case_->fluid_bodies)
                std::cout<<" "<<fb.n<<" = "<<case_->m.body(fb.body).name;
            std::cout<<"), loads of order "<<case_->fluid_order<<(case_->fluid_predictor ? ", predictor" : "")<<std::endl;
            if(case_->run.static_ini)
                std::cout<<"CORAL: static equilibrium with fluid bodies: the fluid loads are not known yet (zero)"<<std::endl;
        }
    }
}

coral_reef3d::~coral_reef3d()
{
    if(sixdof_external::instance()==this)
        sixdof_external::attach(nullptr);
    if(out_)
        out_->close();
}

void coral_reef3d::fail(ghostcell *pgc, const std::string &msg)
{
    std::cerr<<"\n!!! REEF3D::CORAL: "<<msg<<" !!!\n"<<std::endl;
    MPI_Abort(pgc->mpi_comm, 1);
}

void coral_reef3d::ini(lexer *p, ghostcell*)
{
    t_ = p->simtime;
    if(case_->run.static_ini)
    {
        const coral::static_result sr = coral::static_solve(case_->m, t_);
        if(p->mpirank==0)
            std::cout<<"CORAL static equilibrium: "<<(sr.converged ? "converged" : "NOT converged")<<" after "
                     <<sr.iterations<<" iterations"<<std::endl;
    }
    if(p->mpirank==0)
        mkdir("./REEF3D_CORAL", 0777);
    if(case_->run.modes > 0)
    {
        const coral::mode_result mr = coral::modes(case_->m, t_, case_->run.modes, case_->run.modes_damped);
        if(p->mpirank==0)
            coral::write_modes("./REEF3D_CORAL/", case_->m, mr);
    }
    ti_->initialize(case_->m, t_);
    if(p->mpirank==0)
    {
        std::string f = case_->run.out_file;
        const size_t s = f.find_last_of('/');
        if(s != std::string::npos)
            f = f.substr(s+1);
        out_ = std::make_unique<coral::time_series>();
        out_->open("./REEF3D_CORAL/" + f, case_->m);
        out_->write(case_->m, t_);
    }
    next_out_ = t_ + case_->run.out_dt;
}

void coral_reef3d::step(lexer *p, ghostcell*)
{
    const double t1 = p->simtime + p->dt;
    // the fluid bodies have already advanced CORAL within the fluid step
    if(t1 <= t_ + 1.0e-12*std::max(1.0, std::fabs(t1)))
        return;
    if(case_->fluid && !fluid_started_ && nstep_ == 0 && p->mpirank == 0)
        std::cout<<"CORAL: WARNING fluid_body in coral.txt, but no 6DOF body of the fluid solver is connected "
                   "(FNPF with X 10 1 so far); the fluid loads are zero"<<std::endl;
    advance_to(t1);
}

void coral_reef3d::advance_to(double t1)
{
    const double hmax = case_->run.dt > 0.0 ? case_->run.dt : t1 - t_;
    const int n = std::max(1, (int)std::ceil((t1 - t_)/hmax - 1.0e-9));
    const double h = (t1 - t_)/n;
    track_.clear();
    track_.record(case_->m, t_);
    for(int k=0; k<n; ++k)
    {
        const double t = t_ + k*h;
        if(!ti_->step(case_->m, t, h))
        {
            // sub-steps of the failed step
            bool ok = false;
            for(int level=1; level<=6 && !ok; ++level)
            {
                const int m = 1<<level;
                ok = true;
                for(int s=0; s<m && ok; ++s)
                    ok = ti_->step(case_->m, t + s*h/m, h/m);
            }
            if(!ok)
                fail(pgc_, "time step failed at t = " + std::to_string(t) + " s");
        }
        ++nstep_;
        track_.record(case_->m, t + h);
        if(out_ && (case_->run.out_dt <= 0.0 || t + h >= next_out_ - 1.0e-9*h))
        {
            out_->write(case_->m, t + h);
            if(case_->run.out_dt > 0.0)
                while(next_out_ <= t + h + 1.0e-9*h)
                    next_out_ += case_->run.out_dt;
        }
    }
    t_ = t1;
}

// ---- external dynamics of the 6DOF bodies (6DOF_external.h)

int coral_reef3d::fluid_index(int n) const
{
    for(size_t k=0; k<case_->fluid_bodies.size(); ++k)
        if(case_->fluid_bodies[k].n == n)
            return (int)k;
    return -1;
}

bool coral_reef3d::drives(int n) const
{
    return fluid_index(n) >= 0;
}

bool coral_reef3d::predictor() const
{
    return case_->fluid_predictor;
}

void coral_reef3d::connect(int n, double mass, const Eigen::Vector3d &c, const Eigen::Matrix3d &R, const Eigen::Matrix3d &I)
{
    const int k = fluid_index(n);
    if(k < 0)
        return;
    const auto &spec = case_->fluid_bodies[k];
    coral::rigid_body &b = case_->m.bodies[spec.body];
    // mass and inertia of the REEF3D body (its inertia about the CoG, rotated into CORAL's body frame)
    if(!spec.own_mass)
    {
        b.mass = mass;
        b.inertia = b.R.transpose()*(R*I*R.transpose())*b.R;
    }
    const double dx = (b.x - c).norm();
    const double dr = coral::log_so3(b.R*R.transpose()).norm();
    if(p_->mpirank==0)
    {
        std::cout<<"CORAL: fluid body "<<b.name<<" = 6DOF body "<<n<<": mass "<<b.mass<<" kg"
                 <<(spec.own_mass ? " (coral.txt)" : " (REEF3D)")<<", CoG CORAL "<<b.x.transpose()<<" / REEF3D "<<c.transpose()<<std::endl;
        if(dx > 1.0e-3 || dr > 1.0e-3)
            std::cout<<"CORAL: WARNING fluid body "<<b.name<<": CORAL and REEF3D start "<<dx<<" m and "<<dr
                     <<" rad apart, CORAL's pose is used (cog and rpy in coral.txt)"<<std::endl;
        const coral::Vec3 g(p_->W20, p_->W21, p_->W22);
        if((case_->m.gravity - g).norm() > 1.0e-6*std::max(1.0, g.norm()))
            std::cout<<"CORAL: WARNING gravity "<<case_->m.gravity.transpose()<<" in coral.txt, "<<g.transpose()
                     <<" in ctrl.txt (W 20-22); the fluid loads come without the weight"<<std::endl;
        if(spec.own_mass && std::fabs(b.mass - mass) > 1.0e-3*mass)
            std::cout<<"CORAL: WARNING fluid body "<<b.name<<": mass "<<b.mass<<" kg in coral.txt, "<<mass<<" kg in REEF3D"<<std::endl;
    }
}

void coral_reef3d::fluid_loads(double t, const std::vector<vec6> &F, const std::vector<mat6> &A, bool provisional)
{
    const int N = (int)case_->fluid_bodies.size();
    coral::VecX f = coral::VecX::Zero(6*N);
    coral::MatX Am;
    for(int k=0; k<N; ++k)
    {
        const int n = case_->fluid_bodies[k].n;
        if(n < (int)F.size())
            f.segment<6>(6*k) = F[n];
    }
    const int nr = (int)std::lround(std::sqrt((double)A.size()));
    if(!A.empty() && nr*nr == (int)A.size())
    {
        Am = coral::MatX::Zero(6*N, 6*N);
        for(int k=0; k<N; ++k)
            for(int l=0; l<N; ++l)
            {
                const int nk = case_->fluid_bodies[k].n, nl = case_->fluid_bodies[l].n;
                if(nk < nr && nl < nr)
                    Am.block<6,6>(6*k, 6*l) = A[nk*nr + nl];
            }
    }
    case_->fluid->set_sample(t, f, Am, provisional);
    // first loads: consistent accelerations of the start with the fluid loads, masses of the fluid
    if(!fluid_started_ && !provisional)
    {
        ti_->initialize(case_->m, t_);
        fluid_started_ = true;
    }
}

void coral_reef3d::advance(double t0, double t1)
{
    if(std::fabs(t0 - t_) > 1.0e-9*std::max(1.0, std::fabs(t0)) && p_->mpirank==0)
        std::cout<<"CORAL: WARNING fluid step from "<<t0<<" s, CORAL is at "<<t_<<" s"<<std::endl;
    if(t1 > t_)
        advance_to(t1);
}

void coral_reef3d::predict(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w)
{
    // Taylor series from the current state (position second order, rotation with the angular velocity
    // and acceleration, velocities first order)
    const int k = fluid_index(n);
    const coral::rigid_body &b = case_->m.body(case_->fluid_bodies[std::max(0, k)].body);
    const double tau = t - t_;
    const coral::Vec3 wg = b.R*b.w, ag = b.R*b.wd;
    c = b.x + tau*b.v + 0.5*tau*tau*b.vd;
    v = b.v + tau*b.vd;
    R = coral::exp_so3(tau*wg + 0.5*tau*tau*ag)*b.R;
    w = wg + tau*ag;
}

void coral_reef3d::state(int n, double t, Eigen::Vector3d &c, Eigen::Matrix3d &R, Eigen::Vector3d &v, Eigen::Vector3d &w)
{
    const int k = std::max(0, fluid_index(n));
    if(track_.size() == 0)
    {
        predict(n, t, c, R, v, w);
        return;
    }
    const coral::body_snapshot q = track_.state(k, t);
    c = q.x;
    R = q.R;
    v = q.v;
    w = q.w;
}
