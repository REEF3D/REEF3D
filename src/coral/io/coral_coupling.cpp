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

#include"coral/io/coral_coupling.h"
#include"coral/core/coral_state.h"
#include"coral/core/coral_static.h"
#include"coral/core/coral_rotation.h"
#include"coral/linear/coral_modes.h"
#include"coral/io/coral_modes_out.h"
#include"coral/hydro/coral_fluid_load.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include"coral/hydro/coral_member_fluid_load.h"
#include<algorithm>
#include<cmath>
#include<cstdint>
#include<stdexcept>

namespace coral
{

namespace
{
    // journal events: tag, count, doubles
    enum jtag {J_KIN = 1, J_INI = 2, J_CONNECT = 3, J_LOADS = 4, J_STRIPS = 5, J_STEPTO = 6, J_ADVANCE = 7, J_KIN_BED = 8};
    const char J_MAGIC[8] = {'C','O','R','A','L','J','1','\n'};

    std::string file_name(const std::string &f)
    {
        const size_t s = f.find_last_of('/');
        return s == std::string::npos ? f : f.substr(s+1);
    }
}

coupling::coupling(coral_case &c, std::ostream &log, bool out)
    : c_(c), log_(log), out_on_(out)
{
    ti_ = std::make_unique<timeint>(c_.run.ti);
    std::vector<int> b;
    for(const auto &fb : c_.fluid_bodies)
    {
        b.push_back(fb.body);
        if(!fb.one_way)
            twoway_ = true;
    }
    track_.set_bodies(b);
}

coupling::~coupling()
{
    if(out_)
        out_->close();
}

void coupling::jrec(int tag, const std::vector<double> &v)
{
    if(!jout_.is_open())
        return;
    const int32_t t = tag;
    const int64_t n = (int64_t)v.size();
    jout_.write(reinterpret_cast<const char*>(&t), sizeof(t));
    jout_.write(reinterpret_cast<const char*>(&n), sizeof(n));
    if(n)
        jout_.write(reinterpret_cast<const char*>(v.data()), n*sizeof(double));
}

void coupling::record(const std::string &dir)
{
    if(!out_on_ || c_.fluid_record.empty())
        return;
    const std::string f = dir + "/" + file_name(c_.fluid_record);
    // a restart continues the journal (the events from the checkpoint on)
    if(restart_)
    {
        std::ifstream test(f, std::ios::binary);
        if(test.good())
        {
            jout_.open(f, std::ios::binary | std::ios::app);
            if(!jout_)
                throw std::runtime_error("cannot write the journal " + f);
            log_<<"CORAL: journal of the fluid input continued in "<<f<<std::endl;
            return;
        }
    }
    jout_.open(f, std::ios::binary);
    if(!jout_)
        throw std::runtime_error("cannot write the journal " + f);
    jout_.write(J_MAGIC, 8);
    log_<<"CORAL: journal of the fluid input in "<<f<<std::endl;
}

void coupling::ini(double t0, const std::string &dir)
{
    if(!restart_)
    jrec(J_INI, {t0});
    t_ = t0;
    // a restart: the state follows (load_state), no static equilibrium and no modes
    if(c_.run.static_ini && !restart_)
    {
        const static_result sr = static_solve(c_.m, t_);
        log_<<"CORAL static equilibrium: "<<(sr.converged ? "converged" : "NOT converged")<<" after "<<sr.iterations<<" iterations"<<std::endl;
    }
    if(c_.run.modes > 0 && !restart_)
    {
        const mode_result mr = modes(c_.m, t_, c_.run.modes, c_.run.modes_damped);
        if(out_on_)
            write_modes(dir + "/", c_.m, mr);
    }
    ti_->initialize(c_.m, t_);
    if(out_on_)
    {
        out_ = std::make_unique<time_series>();
        out_->open(dir + "/" + file_name(c_.run.out_file), c_.m, restart_);
        if(!restart_)
        out_->write(c_.m, t_);
    }
    next_out_ = t_ + c_.run.out_dt;
}

double coupling::step_to(double t1)
{
    // the fluid bodies have already advanced CORAL within the fluid step
    if(t1 <= t_ + 1.0e-12*std::max(1.0, std::fabs(t1)))
        return t_;
    // one-way loads only: CORAL follows the newest loads (FNPF: those of the start of the step, so CORAL
    // trails the fluid by one step; NHFLOW: those of its end)
    double target = t1;
    const bool fluid = c_.fluid || !c_.fluid_members.empty();
    if(fluid && !twoway_ && fluid_started_)
    {
        target = std::min(t1, t_load_);
        if(target <= t_ + 1.0e-12*std::max(1.0, std::fabs(t1)))
            return t_;
    }
    if(fluid && !fluid_started_ && nstep_ == 0)
        log_<<"CORAL: WARNING fluid_body or fluid_member in coral.txt, but no 6DOF body of the fluid solver is connected "
              "(FNPF, or NHFLOW one-way); the fluid loads are zero"<<std::endl;
    jrec(J_STEPTO, {target});
    advance_to(target);
    return t_;
}

void coupling::advance(double t0, double t1)
{
    jrec(J_ADVANCE, {t0, t1});
    if(std::fabs(t0 - t_) > 1.0e-9*std::max(1.0, std::fabs(t0)))
        log_<<"CORAL: WARNING fluid step from "<<t0<<" s, CORAL is at "<<t_<<" s"<<std::endl;
    if(t1 > t_)
        advance_to(t1);
}

void coupling::advance_to(double t1)
{
    const double hmax = c_.run.dt > 0.0 ? c_.run.dt : t1 - t_;
    const int n = std::max(1, (int)std::ceil((t1 - t_)/hmax - 1.0e-9));
    const double h = (t1 - t_)/n;
    track_.clear();
    track_.record(c_.m, t_);
    for(int k=0; k<n; ++k)
    {
        const double t = t_ + k*h;
        if(!ti_->step(c_.m, t, h))
        {
            // sub-steps of the failed step
            bool ok = false;
            for(int level=1; level<=6 && !ok; ++level)
            {
                const int m = 1<<level;
                ok = true;
                for(int s=0; s<m && ok; ++s)
                    ok = ti_->step(c_.m, t + s*h/m, h/m);
            }
            if(!ok)
                throw std::runtime_error("time step failed at t = " + std::to_string(t) + " s");
        }
        ++nstep_;
        track_.record(c_.m, t + h);
        if(out_ && (c_.run.out_dt <= 0.0 || t + h >= next_out_ - 1.0e-9*h))
        {
            out_->write(c_.m, t + h);
            if(c_.run.out_dt > 0.0)
                while(next_out_ <= t + h + 1.0e-9*h)
                    next_out_ += c_.run.out_dt;
        }
    }
    t_ = t1;
    if(jout_.is_open())
        jout_.flush();
}

// ---- the bodies of the fluid solver

int coupling::fluid_index(int n) const
{
    for(size_t k=0; k<c_.fluid_bodies.size(); ++k)
        if(c_.fluid_bodies[k].n == n)
            return (int)k;
    return -1;
}

int coupling::member_index(int n) const
{
    for(size_t k=0; k<c_.fluid_members.size(); ++k)
        if(c_.fluid_members[k].n == n)
            return (int)k;
    return -1;
}

bool coupling::drives(int n) const
{
    return fluid_index(n) >= 0 || member_index(n) >= 0;
}

bool coupling::one_way(int n) const
{
    const int k = fluid_index(n);
    return k >= 0 ? c_.fluid_bodies[k].one_way : member_index(n) >= 0;
}

bool coupling::strips(int n, Vec3 &a, Vec3 &e, double &L, int &ns) const
{
    const int k = member_index(n);
    if(k < 0)
        return false;
    const member_fluid_load *l = c_.fluid_members[k].load;
    a = l->p1() + c_.fluid_origin;
    e = l->axis();
    L = l->length();
    ns = l->nstrips();
    return true;
}

void coupling::connect(int n, double mass, const Vec3 &c, const Mat3 &R, const Mat3 &I)
{
    // a restart from a checkpoint: connected in the run that wrote it
    if(restart_)
        return;
    {
        std::vector<double> v = {(double)n, mass, c(0), c(1), c(2)};
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) v.push_back(R(i,j));
        for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) v.push_back(I(i,j));
        jrec(J_CONNECT, v);
    }
    const Vec3 o = c_.fluid_origin;
    const int km = member_index(n);
    if(km >= 0)
    {
        const member_fluid_load *l = c_.fluid_members[km].load;
        log_<<"CORAL: loads of 6DOF body "<<n<<" (CoG "<<c.transpose()<<", fluid frame) in "<<l->nstrips()
            <<" strips along "<<l->name()<<", "<<l->length()<<" m from "<<(l->p1() + o).transpose()<<std::endl;
    }
    const int k = fluid_index(n);
    if(k < 0)
        return;
    const auto &spec = c_.fluid_bodies[k];
    rigid_body &b = c_.m.bodies[spec.body];
    // one-way: the loads act at the REEF3D CoG on the CORAL body, nothing else is taken over
    if(spec.one_way)
    {
        const Vec3 pc = c - o;
        c_.fluid->set_attachment(c_.m, k, pc);
        log_<<"CORAL: one-way loads of 6DOF body "<<n<<" on "<<b.name<<" at "<<pc.transpose()
            <<(spec.excitation ? " (without the loads of the start)" : "")<<std::endl;
        return;
    }
    // mass and inertia of the REEF3D body (its inertia about the CoG, rotated into CORAL's body frame)
    if(!spec.own_mass)
    {
        b.mass = mass;
        b.inertia = b.R.transpose()*(R*I*R.transpose())*b.R;
    }
    const double dx = (b.x + o - c).norm();
    const double dr = log_so3(b.R*R.transpose()).norm();
    log_<<"CORAL: fluid body "<<b.name<<" = 6DOF body "<<n<<": mass "<<b.mass<<" kg"
        <<(spec.own_mass ? " (coral.txt)" : " (REEF3D)")<<", CoG CORAL "<<(b.x + o).transpose()<<" / REEF3D "<<c.transpose()<<" (fluid frame)"<<std::endl;
    if(dx > 1.0e-3 || dr > 1.0e-3)
        log_<<"CORAL: WARNING fluid body "<<b.name<<": CORAL and REEF3D start "<<dx<<" m and "<<dr
            <<" rad apart, CORAL's pose is used (cog and rpy in coral.txt)"<<std::endl;
    if(spec.own_mass && std::fabs(b.mass - mass) > 1.0e-3*mass)
        log_<<"CORAL: WARNING fluid body "<<b.name<<": mass "<<b.mass<<" kg in coral.txt, "<<mass<<" kg in REEF3D"<<std::endl;
}

void coupling::fluid_strips(int n, double t, const VecX &F)
{
    {
        std::vector<double> v = {(double)n, t};
        v.insert(v.end(), F.data(), F.data() + F.size());
        jrec(J_STRIPS, v);
    }
    const int k = member_index(n);
    if(k >= 0)
        c_.fluid_members[k].load->set_sample(t, F);
}

void coupling::fluid_loads(double t, const std::vector<Vec6> &F, const std::vector<Mat6> &A, bool provisional)
{
    {
        std::vector<double> v = {t, provisional ? 1.0 : 0.0, (double)F.size()};
        for(const auto &f : F) v.insert(v.end(), f.data(), f.data() + 6);
        v.push_back((double)A.size());
        for(const auto &a : A) for(int i=0; i<6; ++i) for(int j=0; j<6; ++j) v.push_back(a(i,j));
        jrec(J_LOADS, v);
    }
    if(!provisional)
        t_load_ = t;
    if(c_.fluid)
    {
        const int N = (int)c_.fluid_bodies.size();
        VecX f = VecX::Zero(6*N);
        MatX Am;
        for(int k=0; k<N; ++k)
        {
            const int n = c_.fluid_bodies[k].n;
            if(n < (int)F.size())
                f.segment<6>(6*k) = F[n];
        }
        const int nr = (int)std::lround(std::sqrt((double)A.size()));
        if(!A.empty() && nr*nr == (int)A.size())
        {
            Am = MatX::Zero(6*N, 6*N);
            for(int k=0; k<N; ++k)
                for(int l=0; l<N; ++l)
                {
                    const int nk = c_.fluid_bodies[k].n, nl = c_.fluid_bodies[l].n;
                    if(nk < nr && nl < nr)
                        Am.block<6,6>(6*k, 6*l) = A[nk*nr + nl];
                }
        }
        // one-way bodies: no added mass of their own accelerations (the fluid does not see them move);
        // excitation: without the loads of the start
        for(int k=0; k<N; ++k)
        {
            const auto &spec = c_.fluid_bodies[k];
            if(spec.one_way && Am.size())
                Am.middleCols<6>(6*k).setZero();
            if(spec.excitation && !provisional)
            {
                if((int)F_ref_.size() != N)
                    F_ref_.assign(N, VecX());
                if(F_ref_[k].size() == 0)
                    F_ref_[k] = f.segment<6>(6*k);
            }
            if(spec.excitation && k < (int)F_ref_.size() && F_ref_[k].size())
                f.segment<6>(6*k) -= F_ref_[k];
        }
        // the accelerations of the bodies the fluid saw in the step of this sample (artificial added mass):
        // committed: the mean of CORAL's motion over the step (v(t) - v(t_last))/(t - t_last), CORAL has
        // stepped to t; provisional: the predicted pose keeps the acceleration of the start
        VecX acc;
        if(c_.fluid->has_artificial_added_mass())
        {
            acc = c_.fluid->acceleration(c_.m);
            if(!provisional)
            {
                const VecX v = c_.fluid->velocity(c_.m);
                if(v_commit_.size() == v.size() && t > t_commit_)
                    acc = (v - v_commit_)/(t - t_commit_);
                v_commit_ = v;
                t_commit_ = t;
            }
        }
        else if(!am_told_ && Am.size() == 0)
            for(const auto &spec : c_.fluid_bodies)
                if(!spec.one_way)
                {
                    log_<<"CORAL: the fluid gives no added mass and fluid_body has no added_mass: the coupling is"
                        <<" explicit in the loads, light bodies may need added_mass (coral/hydro/coral_fluid_load.h)"<<std::endl;
                    break;
                }
        am_told_ = true;
        c_.fluid->set_sample(t, f, Am, provisional, acc);
    }
    // first loads: consistent accelerations of the start with the fluid loads, masses of the fluid
    if(!fluid_started_ && !provisional)
    {
        ti_->initialize(c_.m, t_);
        fluid_started_ = true;
    }
}

void coupling::predict(int n, double t, Vec3 &c, Mat3 &R, Vec3 &v, Vec3 &w) const
{
    // Taylor series from the current state (position second order, rotation with the angular velocity
    // and acceleration, velocities first order)
    const int k = fluid_index(n);
    const rigid_body &b = c_.m.body(c_.fluid_bodies[std::max(0, k)].body);
    const double tau = t - t_;
    const Vec3 wg = b.R*b.w, ag = b.R*b.wd;
    c = b.x + c_.fluid_origin + tau*b.v + 0.5*tau*tau*b.vd;
    v = b.v + tau*b.vd;
    R = exp_so3(tau*wg + 0.5*tau*tau*ag)*b.R;
    w = wg + tau*ag;
}

void coupling::state(int n, double t, Vec3 &c, Mat3 &R, Vec3 &v, Vec3 &w) const
{
    const int k = std::max(0, fluid_index(n));
    if(track_.size() == 0)
    {
        predict(n, t, c, R, v, w);
        return;
    }
    const body_snapshot q = track_.state(k, t);
    c = q.x + c_.fluid_origin;
    R = q.R;
    v = q.v;
    w = q.w;
}

// ---- sampled kinematics

void coupling::kinematics_points(std::vector<Vec3> &pts) const
{
    pts.clear();
    if(!c_.fluid_kin)
        return;
    c_.fluid_kin->collect(c_.m, pts);
    for(auto &p : pts)
        p += c_.fluid_origin;
}

void coupling::kinematics(double t, const std::vector<Vec3> &u, const std::vector<double> &eta, const std::vector<int> &in,
                          const std::vector<double> &bed)
{
    if(!c_.fluid_kin)
        return;
    const size_t n = u.size();
    const bool hb = bed.size()==n && n > 0;
    {
        std::vector<double> v = {t, (double)n};
        for(size_t q=0; q<n; ++q)
        {
            v.push_back(u[q](0)); v.push_back(u[q](1)); v.push_back(u[q](2));
            v.push_back(eta[q]);
            v.push_back(in[q]);
            if(hb) v.push_back(bed[q]);
        }
        jrec(hb ? J_KIN_BED : J_KIN, v);
    }
    std::vector<double> e(n), b;
    for(size_t q=0; q<n; ++q)
        e[q] = eta[q] - c_.fluid_origin(2);
    if(hb)
    {
        b.resize(n);
        for(size_t q=0; q<n; ++q)
            b[q] = bed[q] - c_.fluid_origin(2);
    }
    c_.fluid_kin->add_sample(t, u, e, in, b);
}

// ---- replay

long coupling::replay(coral_case &c, const std::string &journal, const std::string &dir, std::ostream &log,
                      const std::function<void(double)> &after)
{
    std::ifstream in(journal, std::ios::binary);
    char magic[8];
    if(!in || !in.read(magic, 8) || !std::equal(magic, magic + 8, J_MAGIC))
        throw std::runtime_error("'" + journal + "' is not a CORAL journal");
    coupling cp(c, log, true);
    long nev = 0;
    std::vector<double> v;
    while(true)
    {
        int32_t tag;
        int64_t n;
        if(!in.read(reinterpret_cast<char*>(&tag), sizeof(tag)))
            break;
        if(!in.read(reinterpret_cast<char*>(&n), sizeof(n)) || n < 0)
            throw std::runtime_error("journal '" + journal + "' is cut off");
        v.resize(n);
        if(n && !in.read(reinterpret_cast<char*>(v.data()), n*sizeof(double)))
            throw std::runtime_error("journal '" + journal + "' is cut off");
        ++nev;
        switch(tag)
        {
        case J_INI:
            cp.ini(v.at(0), dir);
            break;
        case J_KIN:
        case J_KIN_BED:
        {
            const int np = (int)v.at(1), w = tag==J_KIN_BED ? 6 : 5;
            std::vector<Vec3> pts, u(np);
            std::vector<double> eta(np), bed;
            std::vector<int> ins(np);
            cp.kinematics_points(pts);
            for(int q=0; q<np; ++q)
            {
                u[q] = Vec3(v.at(2+w*q), v.at(3+w*q), v.at(4+w*q));
                eta[q] = v.at(5+w*q);
                ins[q] = (int)v.at(6+w*q);
                if(w==6) bed.push_back(v.at(7+w*q));
            }
            cp.kinematics(v[0], u, eta, ins, bed);
            break;
        }
        case J_CONNECT:
        {
            Mat3 R, I;
            for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) {R(i,j) = v.at(5+3*i+j); I(i,j) = v.at(14+3*i+j);}
            cp.connect((int)v.at(0), v.at(1), Vec3(v.at(2), v.at(3), v.at(4)), R, I);
            break;
        }
        case J_LOADS:
        {
            const int nf = (int)v.at(2);
            std::vector<Vec6> F(nf);
            for(int k=0; k<nf; ++k) for(int i=0; i<6; ++i) F[k](i) = v.at(3+6*k+i);
            const int na = (int)v.at(3+6*nf);
            std::vector<Mat6> A(na);
            for(int k=0; k<na; ++k) for(int i=0; i<6; ++i) for(int j=0; j<6; ++j) A[k](i,j) = v.at(4+6*nf+36*k+6*i+j);
            cp.fluid_loads(v[0], F, A, v.at(1) > 0.5);
            break;
        }
        case J_STRIPS:
        {
            VecX F((int)n - 2);
            for(int i=0; i<F.size(); ++i) F(i) = v[2+i];
            cp.fluid_strips((int)v.at(0), v.at(1), F);
            break;
        }
        case J_STEPTO:
        {
            const long n0 = cp.steps();
            cp.advance_to(v.at(0));
            if(after && cp.steps() > n0)
                after(cp.time());
            break;
        }
        case J_ADVANCE:
        {
            const long n0 = cp.steps();
            cp.advance(v.at(0), v.at(1));
            if(after && cp.steps() > n0)
                after(cp.time());
            break;
        }
        default:
            throw std::runtime_error("journal '" + journal + "': unknown event " + std::to_string(tag));
        }
    }
    return nev;
}

// ---- restart from a checkpoint (coral_state.h)

void coupling::flush()
{
    if(out_)
        out_->flush();
    if(jout_.is_open())
        jout_.flush();
}

void coupling::save_state(std::vector<double> &v) const
{
    state_out o{v};
    o.d(1.0);       // version
    for(double x : {t_, next_out_, t_commit_, t_load_})
        o.d(x);
    o.i(nstep_);
    o.b(am_told_);
    o.b(fluid_started_);
    o.m(v_commit_);
    o.list(F_ref_, [&](const VecX &x){o.m(x);});
    track_.save_state(v);
    ti_->save_state(v);
    c_.m.save_state(v);
    o.b(c_.fluid_kin != nullptr);
    if(c_.fluid_kin)
        c_.fluid_kin->save_state(v);
}

bool coupling::load_state(const std::vector<double> &v)
{
    state_in in{v};
    double version = 0.0;
    in.d(version);
    if(!in.ok || int(version) != 1)
        return false;
    in.d(t_); in.d(next_out_); in.d(t_commit_); in.d(t_load_);
    in.i(nstep_);
    in.b(am_told_);
    in.b(fluid_started_);
    in.m(v_commit_);
    in.list(F_ref_, [&](VecX &x){in.m(x);});
    if(!in.ok)
        return false;
    size_t k = in.k;
    if(!track_.load_state(v, k) || !ti_->load_state(c_.m, v, k) || !c_.m.load_state(v, k))
        return false;
    state_in in2{v, k};
    bool kin = false;
    in2.b(kin);
    if(!in2.ok || kin != (c_.fluid_kin != nullptr))
        return false;
    k = in2.k;
    if(kin && !c_.fluid_kin->load_state(v, k))
        return false;
    return true;
}

}
