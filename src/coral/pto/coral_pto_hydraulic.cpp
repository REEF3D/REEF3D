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

#include"coral/pto/coral_pto_hydraulic.h"
#include<Eigen/Dense>
#include<algorithm>
#include<cmath>

namespace coral
{

std::string hydraulic_settings::check() const
{
    if(A_a <= 0.0 || A_b <= 0.0) return "the piston areas A_a and A_b must be positive";
    if(V_a0 <= 0.0 || V_b0 <= 0.0) return "the chamber volumes V_a0 and V_b0 must be positive";
    if(beta <= 0.0 || rho <= 0.0) return "beta and rho must be positive";
    if(Cd <= 0.0 || A_v <= 0.0) return "the check valves need Cd > 0 and A_v > 0";
    if(dp_open <= 0.0 || dp_crack < 0.0) return "dp_open must be positive, dp_crack >= 0";
    if(V_hp <= 0.0 || V_lp <= 0.0) return "the accumulator volumes V_hp and V_lp must be positive";
    if(p_hp_pre <= 0.0 || p_lp_pre <= 0.0) return "the precharge pressures p_hp_pre and p_lp_pre must be positive";
    if(p_hp0 < p_hp_pre || p_lp0 < p_lp_pre) return "the initial pressures p_hp0, p_lp0 must not be below the precharge";
    if(kappa < 1.0) return "kappa must be >= 1";
    if(D <= 0.0) return "the motor displacement D must be positive";
    if(eta_m <= 0.0 || eta_m > 1.0) return "eta_m must be in (0, 1]";
    if(J <= 0.0 || c_w < 0.0 || C_leak < 0.0 || b_gen < 0.0 || K_r < 0.0) return "J > 0 and c_w, C_leak, b_gen, K_r >= 0 needed";
    if(dt_sub <= 0.0) return "dt_sub must be positive";
    if(has_gen)
    {
        const std::string g = generator_model(gen).check();
        if(!g.empty()) return "generator: " + g;
    }
    return "";
}

pto_hydraulic::pto_hydraulic(const hydraulic_settings &s) : set(s), gen_(s.gen)
{
    // the initial state (also before start(): an element used without it starts at its first call)
    yc_[0] = yc_[1] = s.p_lp0;
    yc_[2] = s.V_hp*std::pow(s.p_hp_pre/s.p_hp0, 1.0/s.kappa);
    yc_[3] = s.V_lp*std::pow(s.p_lp_pre/s.p_lp0, 1.0/s.kappa);
    yc_[4] = s.w0;
    y_ = yc_;
}

double pto_hydraulic::p_hp(double Vg) const
{
    return set.p_hp_pre*std::pow(set.V_hp/std::max(Vg, 1.0e-9*set.V_hp), set.kappa);
}

double pto_hydraulic::p_lp(double Vg) const
{
    return set.p_lp_pre*std::pow(set.V_lp/std::max(Vg, 1.0e-9*set.V_lp), set.kappa);
}

double pto_hydraulic::valve(double dp) const
{
    // check valve: flow in its direction for dp > dp_crack, smooth opening over dp_open
    if(dp <= set.dp_crack)
        return 0.0;
    const double u = std::min(1.0, (dp - set.dp_crack)/set.dp_open);
    const double open = u*u*(3.0 - 2.0*u);
    return set.Cd*set.A_v*open*std::sqrt(2.0*dp/set.rho);
}

double pto_hydraulic::gen_torque(double w) const
{
    double T = set.b_gen*w;
    if(set.T_max > 0.0)
        T = std::clamp(T, -set.T_max, set.T_max);
    return T;
}

pto_hydraulic::flows pto_hydraulic::eval_flows(double, const state &y) const
{
    flows f;
    const double ph = p_hp(y[2]), pl = p_lp(y[3]);
    f.Qa_hp = valve(y[0] - ph);
    f.Qb_hp = valve(y[1] - ph);
    f.Qlp_a = valve(pl - y[0]);
    f.Qlp_b = valve(pl - y[1]);
    const double dp = ph - pl, w = y[4];
    f.Qm = set.D*w + set.C_leak*dp;
    f.Qr = (set.K_r > 0.0) ? set.K_r*std::max(0.0, ph - set.p_relief) : 0.0;
    f.Tm = (dp*w >= 0.0) ? set.eta_m*set.D*dp : set.D*dp/set.eta_m;
    f.Tg = gen_torque(w);
    f.P_diss = (y[0] - ph)*f.Qa_hp + (y[1] - ph)*f.Qb_hp + (pl - y[0])*f.Qlp_a + (pl - y[1])*f.Qlp_b
             + set.C_leak*dp*dp + (set.D*dp - f.Tm)*w + set.c_w*w*w + dp*f.Qr;
    return f;
}

void pto_hydraulic::rhs(double x, double v, const state &y, state &dy, flows *fl) const
{
    const flows f = eval_flows(x, y);
    const double Va = std::max(set.V_a0 + set.A_a*x, 1.0e-3*set.V_a0);
    const double Vb = std::max(set.V_b0 - set.A_b*x, 1.0e-3*set.V_b0);
    dy[0] = set.beta/Va*(-set.A_a*v - f.Qa_hp + f.Qlp_a);
    dy[1] = set.beta/Vb*( set.A_b*v - f.Qb_hp + f.Qlp_b);
    dy[2] = -(f.Qa_hp + f.Qb_hp - f.Qm - f.Qr);
    dy[3] = -(f.Qm + f.Qr - f.Qlp_a - f.Qlp_b);
    dy[4] = (f.Tm - f.Tg - set.c_w*y[4])/set.J;
    if(fl)
        *fl = f;
}

void pto_hydraulic::integrate(const pto_state &s, state &y, double *e) const
{
    y = yc_;
    const double h = s.t - tc_;
    if(!(h > 0.0))
        return;
    const int n = std::max(1, (int)std::ceil(h/set.dt_sub - 1.0e-9));
    const double hs = h/n;
    const double q0 = qc_, v0 = qdc_, q1 = s.q, v1 = s.qd;
    const double pref = std::max({set.p_hp0, set.p_lp0, 1.0e5});
    const std::array<double,5> sc = {pref, pref, set.V_hp, set.V_lp, 1.0};
    // piston on the cubic Hermite path between the accepted and the asked state, u in [0, 1]
    auto path = [&](double u, double &x, double &v)
    {
        const double u2 = u*u, u3 = u2*u;
        x = (2.0*u3 - 3.0*u2 + 1.0)*q0 + (u3 - 2.0*u2 + u)*h*v0 + (-2.0*u3 + 3.0*u2)*q1 + (u3 - u2)*h*v1;
        v = ((6.0*u2 - 6.0*u)*q0 + (3.0*u2 - 4.0*u + 1.0)*h*v0 + (-6.0*u2 + 6.0*u)*q1 + (3.0*u2 - 2.0*u)*h*v1)/h;
    };
    // implicit stage: z - hg f(x, v, z) = r, Newton with a finite-difference Jacobian from the guess z
    Eigen::Matrix<double,5,5> Jm;
    Eigen::Matrix<double,5,1> res, d;
    auto stage = [&](double x, double v, const state &r, double hg, state &z)
    {
        state f, f2;
        for(int it=0; it<50; ++it)
        {
            rhs(x, v, z, f);
            for(int i=0; i<5; ++i)
                res(i) = (z[i] - r[i] - hg*f[i])/sc[i];
            for(int j=0; j<5; ++j)
            {
                state zp = z;
                const double dj = 1.0e-7*std::max(std::fabs(z[j]), sc[j]);
                zp[j] += dj;
                rhs(x, v, zp, f2);
                for(int i=0; i<5; ++i)
                    Jm(i,j) = ((i==j ? 1.0 : 0.0) - hg*(f2[i] - f[i])/dj)*sc[j]/sc[i];
            }
            d = Jm.partialPivLu().solve(-res);
            double dmax = 0.0;
            for(int i=0; i<5; ++i)
            {
                z[i] += d(i)*sc[i];
                dmax = std::max(dmax, std::fabs(d(i)));
            }
            z[2] = std::max(z[2], 1.0e-6*set.V_hp);
            z[3] = std::max(z[3], 1.0e-6*set.V_lp);
            if(dmax < 1.0e-12)
                break;
        }
    };
    // two-stage SDIRK (Alexander), L-stable and second order, stiffly accurate:
    //   Y1 = y + g hs f(t + g hs, Y1),  y+ = Y2 = y + (1 - g) hs f(t + g hs, Y1) + g hs f(t + hs, Y2)
    const double g = 1.0 - 0.5*std::sqrt(2.0);
    for(int k=0; k<n; ++k)
    {
        double x1, w1, x2, w2;
        path((k + g)/n, x1, w1);
        path(double(k + 1)/n, x2, w2);
        const state y0 = y;
        state Y1 = y0;
        stage(x1, w1, y0, g*hs, Y1);
        state f1;
        flows fl1;
        rhs(x1, w1, Y1, f1, &fl1);
        state r2;
        for(int i=0; i<5; ++i)
            r2[i] = y0[i] + (1.0 - g)*hs*f1[i];
        state Y2 = Y1;
        stage(x2, w2, r2, g*hs, Y2);
        y = Y2;
        if(e)
        {
            // powers with the weights of the method
            flows fl2;
            state f2;
            rhs(x2, w2, Y2, f2, &fl2);
            const double b1 = (1.0 - g)*hs, b2 = g*hs;
            e[0] += -piston_force(Y1)*w1*b1 - piston_force(Y2)*w2*b2;
            e[1] += fl1.Tg*Y1[4]*b1 + fl2.Tg*Y2[4]*b2;
            e[2] += fl1.P_diss*b1 + fl2.P_diss*b2;
        }
    }
}

pto_output pto_hydraulic::force(const pto_state &s)
{
    if(!started_)
        start(s);
    if(s.t == t_last_ && s.q == q_last_ && s.qd == qd_last_)
        return out_last_;
    integrate(s, y_, nullptr);
    pto_output o;
    o.F = piston_force(y_);
    o.P = -o.F*s.qd;
    t_last_ = s.t; q_last_ = s.q; qd_last_ = s.qd;
    out_last_ = o;
    return o;
}

void pto_hydraulic::start(const pto_state &s)
{
    started_ = true;
    tc_ = s.t;
    qc_ = s.q;
    qdc_ = s.qd;
    y_ = yc_;
    t_last_ = -1.0e300;
    const generator_power gp = gen_.convert(-gen_torque(yc_[4]), yc_[4]);
    Pel_last_ = set.has_gen ? gp.P_grid : 0.0;
}

void pto_hydraulic::accept(const pto_state &s, double)
{
    if(!started_)
        start(s);
    double e[3] = {0.0, 0.0, 0.0};
    const double dt = s.t - tc_;
    integrate(s, yc_, e);
    E_in += e[0]; E_gen += e[1]; E_diss += e[2];
    tc_ = s.t; qc_ = s.q; qdc_ = s.qd;
    y_ = yc_;
    t_last_ = -1.0e300;
    if(set.has_gen)
    {
        const double P = gen_.convert(-gen_torque(yc_[4]), yc_[4]).P_grid;
        if(dt > 0.0)
            E_el += 0.5*dt*(P + Pel_last_);
        Pel_last_ = P;
    }
}

void pto_hydraulic::set_param(const std::string &key, double val)
{
    if(key == "b_gen")
    {
        set.b_gen = val;
        t_last_ = -1.0e300;
    }
}

double pto_hydraulic::gas_energy(const state &y) const
{
    // work done on the gas from the initial volume: -int p dV
    auto w = [&](double V, double ppre, double V0, double Vg) {
        const double k = set.kappa, c = ppre*std::pow(V, k);
        if(std::fabs(k - 1.0) < 1.0e-12)
            return -c*std::log(Vg/V0);
        return c*(std::pow(Vg, 1.0 - k) - std::pow(V0, 1.0 - k))/(k - 1.0);
    };
    const double Vh0 = set.V_hp*std::pow(set.p_hp_pre/set.p_hp0, 1.0/set.kappa);
    const double Vl0 = set.V_lp*std::pow(set.p_lp_pre/set.p_lp0, 1.0/set.kappa);
    return w(set.V_hp, set.p_hp_pre, Vh0, y[2]) + w(set.V_lp, set.p_lp_pre, Vl0, y[3]);
}

void pto_hydraulic::channel_names(const std::string &p, std::vector<std::string> &c) const
{
    for(const char *x : {"pa", "pb", "php", "plp", "w", "Qm", "Pgen"})
        c.push_back(p + "_" + x);
    if(set.has_gen)
        for(const char *x : {"Pe", "Ee"})
            c.push_back(p + "_" + x);
}

void pto_hydraulic::channel_values(std::vector<double> &v) const
{
    const flows f = eval_flows(0.0, y_);
    v.push_back(y_[0]);
    v.push_back(y_[1]);
    v.push_back(p_hp(y_[2]));
    v.push_back(p_lp(y_[3]));
    v.push_back(y_[4]);
    v.push_back(f.Qm);
    v.push_back(f.Tg*y_[4]);
    if(set.has_gen)
    {
        v.push_back(gen_.convert(-f.Tg, y_[4]).P_grid);
        v.push_back(E_el);
    }
}

void pto_hydraulic::save_state(std::vector<double> &v) const
{
    pto_base::save_state(v);
    for(double x : yc_)
        v.push_back(x);
    for(double x : {double(started_), tc_, qc_, qdc_, Pel_last_, set.b_gen, E_in, E_gen, E_diss, E_el})
        v.push_back(x);
}

void pto_hydraulic::load_state(const std::vector<double> &v, size_t &k)
{
    pto_base::load_state(v, k);
    if(k + 15 > v.size())
    {
        k = v.size();
        return;
    }
    for(double &x : yc_)
        x = v[k++];
    started_ = v[k++] > 0.5;
    tc_ = v[k++]; qc_ = v[k++]; qdc_ = v[k++]; Pel_last_ = v[k++]; set.b_gen = v[k++];
    E_in = v[k++]; E_gen = v[k++]; E_diss = v[k++]; E_el = v[k++];
    y_ = yc_;
    t_last_ = -1.0e300;
}

}
