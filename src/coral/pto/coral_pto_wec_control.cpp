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

#include"coral/pto/coral_pto_wec_control.h"
#include<Eigen/Dense>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>
#include<stdexcept>
#if defined(_WIN32)
#else
#include<dlfcn.h>
#endif

namespace coral
{

// --- matrix exponential ---
MatX coral_expm(const MatX &A)
{
    const double nrm = A.cwiseAbs().rowwise().sum().maxCoeff();
    int s = 0;
    if(nrm > 0.5)
        s = std::max(0, (int)std::ceil(std::log2(nrm/0.5)));
    const MatX B = A/std::ldexp(1.0, s);
    MatX E = MatX::Identity(A.rows(), A.cols()), T = E;
    for(int k=1; k<=18; ++k)
    {
        T = T*B/double(k);
        E += T;
        if(T.cwiseAbs().maxCoeff() < 1.0e-18*E.cwiseAbs().maxCoeff())
            break;
    }
    for(int k=0; k<s; ++k)
        E = E*E;
    return E;
}

// --- gain schedule ---
pto_ctrl_schedule::pto_ctrl_schedule(const std::vector<double> &t, const std::vector<double> &k,
                                     const std::vector<double> &c, const std::vector<double> &F0)
    : t_(t), k_(k), c_(c), F0_(F0)
{
}

double pto_ctrl_schedule::value(const std::vector<double> &f, double t) const
{
    if(t <= t_.front())
        return f.front();
    if(t >= t_.back())
        return f.back();
    const size_t n = size_t(std::upper_bound(t_.begin(), t_.end(), t) - t_.begin());
    const double r = (t - t_[n-1])/(t_[n] - t_[n-1]);
    return (1.0 - r)*f[n-1] + r*f[n];
}

void pto_ctrl_schedule::update(const pto_ctrl_in &in, pto_composite &pto)
{
    pto.set_param("spring-damper", "K", value(k_, in.t));
    pto.set_param("spring-damper", "B", value(c_, in.t));
    pto.set_param("spring-damper", "F0", value(F0_, in.t));
}

bool pto_ctrl_schedule::read(const std::string &file, std::vector<double> &t, std::vector<double> &k,
                             std::vector<double> &c, std::vector<double> &F0, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    t.clear(); k.clear(); c.clear(); F0.clear();
    std::string line;
    int n = 0;
    while(std::getline(in, line))
    {
        ++n;
        const size_t h = line.find('#');
        if(h != std::string::npos) line.erase(h);
        std::istringstream ls(line);
        std::vector<double> v;
        double x;
        while(ls >> x) v.push_back(x);
        if(v.empty())
            continue;
        if(v.size() < 3 || v.size() > 4)
        {
            msg = file + " line " + std::to_string(n) + ": 't k c [F0]' expected";
            return false;
        }
        if(!t.empty() && v[0] <= t.back())
        {
            msg = file + " line " + std::to_string(n) + ": the times must increase";
            return false;
        }
        t.push_back(v[0]); k.push_back(v[1]); c.push_back(v[2]); F0.push_back(v.size() > 3 ? v[3] : 0.0);
    }
    if(t.empty())
    {
        msg = file + ": no line 't k c [F0]'";
        return false;
    }
    return true;
}

// --- user library ---
pto_ctrl_user::pto_ctrl_user(const std::string &lib, const std::vector<double> &params, double period,
                             double K, double B, double F0)
    : lib_(lib), params_(params), period_(period)
{
    out_[0] = K; out_[1] = B; out_[2] = F0;
#if defined(_WIN32)
    throw std::runtime_error("PTO user controller: shared libraries are not supported on this platform");
#else
    handle_ = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL);
    if(!handle_)
    {
        const char *e = dlerror();
        throw std::runtime_error("PTO user controller: cannot load '" + lib + "': " + (e ? e : "unknown error"));
    }
    fn_ = reinterpret_cast<fn_t>(dlsym(handle_, "coral_pto_control"));
    if(!fn_)
    {
        dlclose(handle_);
        handle_ = nullptr;
        throw std::runtime_error("PTO user controller: '" + lib + "' has no function coral_pto_control "
                                 "(declare it extern \"C\" in C++)");
    }
    save_ = reinterpret_cast<save_t>(dlsym(handle_, "coral_pto_control_save"));
    load_ = reinterpret_cast<load_t>(dlsym(handle_, "coral_pto_control_load"));
    end_ = reinterpret_cast<end_t>(dlsym(handle_, "coral_pto_control_end"));
#endif
}

pto_ctrl_user::~pto_ctrl_user()
{
#if !defined(_WIN32)
    if(end_)
        end_(user_);
    if(handle_)
        dlclose(handle_);
#endif
}

void pto_ctrl_user::update(const pto_ctrl_in &in, pto_composite &pto)
{
    if(in.t < t_next_ - 1.0e-9*std::max(1.0, period_))
        return;
    std::vector<double> x = {in.t, in.q, in.qd, in.F, in.P};
    x.insert(x.end(), params_.begin(), params_.end());
    const int rc = fn_((int)x.size(), x.data(), 3, out_, &user_);
    if(rc != 0)
        throw std::runtime_error("PTO user controller '" + lib_ + "' returned " + std::to_string(rc) + " at t = "
                                 + std::to_string(in.t));
    pto.set_param("spring-damper", "K", out_[0]);
    pto.set_param("spring-damper", "B", out_[1]);
    pto.set_param("spring-damper", "F0", out_[2]);
    if(period_ > 0.0)
        t_next_ = (t_next_ < -1.0e299) ? in.t + period_ : t_next_ + period_;
}

void pto_ctrl_user::save_state(std::vector<double> &v) const
{
    for(double x : {t_next_, out_[0], out_[1], out_[2]})
        v.push_back(x);
    int n = 0;
    std::vector<double> b;
    if(save_)
    {
        n = std::max(0, save_(user_, nullptr, 0));
        b.resize(n);
        if(n > 0)
            save_(user_, b.data(), n);
    }
    v.push_back(double(n));
    v.insert(v.end(), b.begin(), b.end());
}

void pto_ctrl_user::load_state(const std::vector<double> &v, size_t &k)
{
    if(k + 5 > v.size())
    {
        k = v.size();
        return;
    }
    t_next_ = v[k++]; out_[0] = v[k++]; out_[1] = v[k++]; out_[2] = v[k++];
    const int n = (int)v[k++];
    if(n > 0 && k + n <= v.size())
    {
        if(load_ && load_(&user_, &v[k], n) != 0)
            throw std::runtime_error("PTO user controller '" + lib_ + "': coral_pto_control_load failed");
        k += n;
    }
}

void pto_ctrl_user::channel_names(const std::string &p, std::vector<std::string> &c) const
{
    for(const char *x : {"uK", "uB", "uF0"})
        c.push_back(p + "_" + x);
}

void pto_ctrl_user::channel_values(std::vector<double> &v) const
{
    v.push_back(out_[0]);
    v.push_back(out_[1]);
    v.push_back(out_[2]);
}

// --- model predictive control ---
pto_ctrl_mpc::pto_ctrl_mpc(const mpc_model &m, const mpc_settings &s) : mod_(m), set_(s)
{
    const int nr = (int)mod_.Ar.rows();
    nx_ = 2 + nr;
    const double Mt = mod_.m + mod_.A_inf;
    if(!(Mt > 0.0))
        throw std::runtime_error("PTO mpc: the mass m + A_inf must be positive");
    if(set_.Tc <= 0.0 || set_.N < 1)
        throw std::runtime_error("PTO mpc: the period Tc must be positive and the horizon N >= 1");
    // continuous model z = [x, v, xi], input u (PTO force) and the excitation enter the same way
    MatX E = MatX::Zero(nx_ + 1, nx_ + 1);
    E(0,1) = 1.0;
    E(1,0) = -mod_.C/Mt;
    E(1,1) = -mod_.B/Mt;
    for(int i=0; i<nr; ++i)
    {
        E(1,2+i) = -mod_.Cr(i)/Mt;
        E(2+i,1) = mod_.Br(i);
        for(int j=0; j<nr; ++j)
            E(2+i,2+j) = mod_.Ar(i,j);
    }
    E(1,nx_) = 1.0/Mt;
    const MatX X = coral_expm(E*set_.Tc);
    Phi_ = X.topLeftCorner(nx_, nx_);
    Gu_ = X.block(0, nx_, nx_, 1);
    xi_ = VecX::Zero(nr);
    u_ = 0.0;
}

void pto_ctrl_mpc::observe(double t, double v)
{
    const int nr = (int)mod_.Ar.rows();
    if(!started_)
    {
        started_ = true;
        t_obs_ = t;
        v_obs_ = v;
        return;
    }
    const double dt = t - t_obs_;
    if(nr > 0 && dt > 0.0)
    {
        if(std::fabs(dt - dt_r_) > 1.0e-12*std::max(1.0, dt))
        {
            MatX E = MatX::Zero(nr + 1, nr + 1);
            E.topLeftCorner(nr, nr) = mod_.Ar;
            E.block(0, nr, nr, 1) = mod_.Br;
            const MatX X = coral_expm(E*dt);
            Phir_ = X.topLeftCorner(nr, nr);
            Gr_ = X.block(0, nr, nr, 1);
            dt_r_ = dt;
        }
        // the velocity linear over the step: its mean held (second order)
        xi_ = Phir_*xi_ + Gr_*(0.5*(v + v_obs_));
    }
    t_obs_ = t;
    v_obs_ = v;
}

VecX pto_ctrl_mpc::solve(double t, double x, double v, const VecX &xi)
{
    const int N = set_.N;
    // free response (u = 0) with the predicted excitation, and the impulse responses of u
    VecX z = VecX::Zero(nx_);
    z(0) = x;
    z(1) = v;
    if(xi.size() == nx_ - 2)
        z.tail(nx_ - 2) = xi;
    VecX av(N+1), ax(N+1);
    av(0) = z(1);
    ax(0) = z(0);
    for(int k=0; k<N; ++k)
    {
        const double Fe = mod_.Fe ? mod_.Fe(t + (k + 0.5)*set_.Tc) : 0.0;
        z = Phi_*z + Gu_*Fe;
        av(k+1) = z(1);
        ax(k+1) = z(0);
    }
    VecX gv(N), gx(N);
    VecX zi = Gu_;
    for(int n=0; n<N; ++n)
    {
        gv(n) = zi(1);
        gx(n) = zi(0);
        zi = Phi_*zi;
    }
    // velocity at the points k = 0..N: v = av + Gv u, Gv(k,j) = gv(k-1-j) for j < k
    MatX Gm = MatX::Zero(N, N), Gx = MatX::Zero(N, N);
    for(int k=0; k<N; ++k)
        for(int j=0; j<=k; ++j)
        {
            const double v0 = (j < k) ? gv(k-1-j) : 0.0;
            Gm(k,j) = 0.5*(v0 + gv(k-j));
            Gx(k,j) = gx(k-j);       // x_{k+1}
        }
    VecX am(N);
    for(int k=0; k<N; ++k)
        am(k) = 0.5*(av(k) + av(k+1));
    MatX H = Gm + Gm.transpose();
    VecX g = am;
    const double dmax = std::max(1.0e-30, H.diagonal().cwiseAbs().maxCoeff());
    // the weight on the force and a regularisation 0.01 Tc/(m + A_inf): the absorbed energy hardly depends on
    // force components at frequencies with little radiation damping (very low and very high), so without it
    // the optimum is not unique there and the forces grow without bound
    H.diagonal().array() += set_.r + 0.01*set_.Tc/(mod_.m + mod_.A_inf);
    // convex: shift by the most negative eigenvalue if there is one
    Eigen::SelfAdjointEigenSolver<MatX> es(H, Eigen::EigenvaluesOnly);
    const double lmin = es.eigenvalues().minCoeff();
    if(lmin < 1.0e-9*dmax)
        H.diagonal().array() += 1.0e-9*dmax - lmin;
    // warm start: the last sequence shifted
    VecX u = VecX::Zero(N);
    if(useq_.size() == N)
    {
        u.head(N-1) = useq_.tail(N-1);
        u(N-1) = useq_(N-1);
    }
    const double Fm = set_.F_max > 0.0 ? set_.F_max : 0.0;
    auto clip = [&](VecX &w) {if(Fm > 0.0) w = w.cwiseMax(-Fm).cwiseMin(Fm);};
    std::vector<int> act;       // stroke penalties
    MatX Hq = H;
    VecX gq = g;
    double wq = set_.w_q;
    iters_ = 0;
    for(int round=0; round<6; ++round)
    {
        // unconstrained minimum; if inside the box it is the solution
        Eigen::LDLT<MatX> ldlt(Hq);
        VecX uu = ldlt.solve(-gq);
        bool inside = (Fm <= 0.0) || uu.cwiseAbs().maxCoeff() <= Fm;
        if(inside)
            u = uu;
        else
        {
            // box: accelerated projected gradient (FISTA) with step 1/L
            Eigen::SelfAdjointEigenSolver<MatX> eh(Hq, Eigen::EigenvaluesOnly);
            const double L = eh.eigenvalues().maxCoeff();
            clip(u);
            VecX y = u, up = u;
            double tk = 1.0;
            for(int it=0; it<5000; ++it)
            {
                VecX un = y - (Hq*y + gq)/L;
                clip(un);
                const double tn = 0.5*(1.0 + std::sqrt(1.0 + 4.0*tk*tk));
                y = un + ((tk - 1.0)/tn)*(un - up);
                const double ch = (un - up).cwiseAbs().maxCoeff();
                up = un;
                tk = tn;
                ++iters_;
                if(ch < 1.0e-10*std::max(Fm, 1.0))
                    break;
            }
            u = up;
        }
        if(set_.q_max <= 0.0)
            break;
        // stroke: quadratic penalty on the predicted excursions beyond q_max
        const VecX xp = ax.tail(N) + Gx*u;
        bool added = false;
        if(wq <= 0.0)
        {
            double gmax = 1.0e-30;
            for(int k=0; k<N; ++k) gmax = std::max(gmax, Gx.row(k).squaredNorm());
            wq = 1.0e3*dmax/gmax;
        }
        for(int k=0; k<N; ++k)
            if(std::fabs(xp(k)) > set_.q_max*(1.0 + 1.0e-6) && std::find(act.begin(), act.end(), k) == act.end())
            {
                act.push_back(k);
                const double sgn = xp(k) > 0.0 ? 1.0 : -1.0;
                Hq += 2.0*wq*Gx.row(k).transpose()*Gx.row(k);
                gq += 2.0*wq*Gx.row(k).transpose()*(ax(k+1) - sgn*set_.q_max);
                added = true;
            }
        if(!added)
            break;
    }
    if(!u.allFinite())
        u.setZero();
    useq_ = u;
    return u;
}

void pto_ctrl_mpc::update(const pto_ctrl_in &in, pto_composite &pto)
{
    observe(in.t, in.qd);
    if(in.t < t_next_ - 1.0e-9*set_.Tc)
        return;
    const VecX u = solve(in.t, in.q, in.qd, xi_);
    u_ = u(0);
    pto.set_param("spring-damper", "K", 0.0);
    pto.set_param("spring-damper", "B", 0.0);
    pto.set_param("spring-damper", "F0", u_);
    t_next_ = (t_next_ < -1.0e299) ? in.t + set_.Tc : t_next_ + set_.Tc;
}

void pto_ctrl_mpc::save_state(std::vector<double> &v) const
{
    for(double x : {double(started_), t_obs_, v_obs_, t_next_, u_, double(xi_.size()), double(useq_.size())})
        v.push_back(x);
    for(int i=0; i<xi_.size(); ++i) v.push_back(xi_(i));
    for(int i=0; i<useq_.size(); ++i) v.push_back(useq_(i));
}

void pto_ctrl_mpc::load_state(const std::vector<double> &v, size_t &k)
{
    if(k + 7 > v.size())
    {
        k = v.size();
        return;
    }
    started_ = v[k++] > 0.5; t_obs_ = v[k++]; v_obs_ = v[k++]; t_next_ = v[k++]; u_ = v[k++];
    const int nr = (int)v[k++], nu = (int)v[k++];
    if(k + nr + nu > v.size())
    {
        k = v.size();
        return;
    }
    xi_.resize(nr);
    for(int i=0; i<nr; ++i) xi_(i) = v[k++];
    useq_.resize(nu);
    for(int i=0; i<nu; ++i) useq_(i) = v[k++];
    dt_r_ = -1.0;
}

void pto_ctrl_mpc::channel_names(const std::string &p, std::vector<std::string> &c) const
{
    c.push_back(p + "_Fmpc");
}

void pto_ctrl_mpc::channel_values(std::vector<double> &v) const
{
    v.push_back(u_);
}

}
