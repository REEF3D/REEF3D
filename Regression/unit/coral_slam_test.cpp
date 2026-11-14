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
// REEF3D::CORAL unit test of the slamming loads of breaking waves (coral/hydro/coral_member_slam, model wo | goda | cw |
// custom): the slamming coefficients at their start, the junction and the end, and their integrals over the impact;
// an impact at a given time on a fixed pile: the peak force rho R V^2 Cs(0) lambda eta_b, its height (the moment about the
// origin), its direction (the crest velocity) and its impulse rho R^2 V lambda eta_b int Cs dtau; the crest trigger with
// the kinematic criterion in a regular wave (one impact per crest, after the crest, at the crest's elevation and
// velocity, none below the threshold); the breaking flag of a sampled fluid (impacts only on the flagged crests); the
// restart state; the keyword and its channels.
// Architect: Hans Bihs
#include"coral/hydro/coral_member_slam.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include"coral/hydro/coral_waves.h"
#include"coral/io/coral_input.h"
#include<cmath>
#include<cstdio>
#include<iostream>
#include<string>
#include<vector>

using namespace coral;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}
static std::string fmt(const char *f, double a, double b=0.0, double c=0.0, double d=0.0)
{
    char s[320];
    std::snprintf(s, sizeof(s), f, a, b, c, d);
    return s;
}
static const double RHO = 1025.0;

// int Cs dtau over the impact (fine midpoint rule)
static double cs_integral(const member_slam_coef &c)
{
    const int n = 400000;
    const double T = slam_duration(c);
    double s = 0.0;
    for(int i=0; i<n; ++i)
        s += slam_coefficient(c, (i + 0.5)*T/n);
    return s*T/n;
}

static void test_coefficients()
{
    std::cout<<"slamming coefficients"<<std::endl;
    member_slam_coef wo, go, cw, cu;
    wo.model = slam_model::wo;
    go.model = slam_model::goda;
    cw.model = slam_model::cw;
    cu.model = slam_model::custom; cu.cs0 = 4.3; cu.tau_end = 0.5;
    const double j1 = slam_coefficient(wo, 0.125), j2 = slam_coefficient(wo, 0.125 + 1.0e-9), we = slam_coefficient(wo, 13.0/32.0);
    check(std::fabs(slam_coefficient(wo, 0.0) - 2.0*M_PI) < 1e-12 && std::fabs(j1 - 4.573) < 2e-3 && std::fabs(j2 - 2.773) < 2e-3
          && std::fabs(we - 1.299) < 2e-3 && slam_coefficient(wo, 13.0/32.0 + 1e-9) == 0.0,
          fmt("wo: 2 pi at the start, %.3f | %.3f at tau = 1/8 (the published form jumps), %.3f at the end 13/32, then 0", j1, j2, we));
    check(std::fabs(slam_coefficient(go, 0.0) - M_PI) < 1e-12 && std::fabs(slam_coefficient(go, 0.5) - 0.5*M_PI) < 1e-12
          && slam_coefficient(go, 1.0 + 1e-9) == 0.0, "goda: pi (1 - tau) until tau = 1");
    check(std::fabs(slam_coefficient(cw, 0.0) - 5.15) < 1e-12 && std::fabs(slam_coefficient(cw, 2.0) - 5.15*(1.0/20.0 + 0.107)) < 1e-12
          && slam_coefficient(cw, 2.0 + 1e-9) == 0.0, "cw: 5.15 (D/(D + 19 s) + 0.107 s/D), 5.15 at the start, until s = D (tau = 2)");
    const double Iw = cs_integral(wo), Ig = cs_integral(go), Ic = cs_integral(cw), Iu = cs_integral(cu);
    check(std::fabs(Ig - 0.5*M_PI) < 1e-6 && std::fabs(Ic - 5.15*(std::log(20.0)/9.5 + 0.107)) < 1e-6 && std::fabs(Iu - 0.5*4.3*0.5) < 1e-6
          && Iw > 0.9 && Iw < 1.3,
          fmt("int Cs dtau: wo %.4f, goda %.4f (pi/2), cw %.4f (5.15 (ln 20/9.5 + 0.107)), custom %.4f", Iw, Ig, Ic, Iu));
}

static member_slam *find_slam(const model &m)
{
    for(const auto &l : m.load_list)
        if(auto *p = dynamic_cast<member_slam*>(l.get()))
            return p;
    return nullptr;
}

// a fixed pile D 6 m from -30 m to +20 m in 1 m elements, in a regular wave
static std::string pile(const std::string &waves, const std::string &slam)
{
    return "rho 1025\ngravity 0 0 -9.81\n" + waves
           + "section pile tube D 6 t 0.06 E 2.1e11 G 8.08e10 rho 7850\n"
             "member mp p1 0 0 -30 p2 0 0 20 n 50 section pile orient 1 0 0 end1 ground\n"
             "slam s member mp " + slam + "\n";
}

static void test_given_impact()
{
    std::cout<<"impact at a given time on a fixed pile"<<std::endl;
    const double V = 10.0, eb = 6.0, lam = 0.46, R = 3.0, t0 = 5.003;
    for(const char *mo : {"wo", "goda", "cw"})
    {
        coral_case c;
        c.read_string(pile("waves regular H 2 T 10 direction 30\n", std::string("model ") + mo + " trigger time at 5.003 eta 6 velocity 10"));
        member_slam *s = find_slam(c.m);
        s->initialize(c.m, 0.0);
        const double dt = 0.01;
        for(int i=1; i<=500; ++i)
            s->commit(c.m, i*dt, dt);
        // the event was created at t = 5.0 for the step to 5.01, starting at 5.003
        check(s->nevents()==1 && std::fabs(s->events()[0].t0 - t0) < 1e-12, std::string(mo) + ": one impact, at its time");
        member_slam_coef cc;
        cc.model = std::string(mo)=="wo" ? slam_model::wo : std::string(mo)=="goda" ? slam_model::goda : slam_model::cw;
        const double Fpk = RHO*R*V*V*slam_coefficient(cc, 0.0)*lam*eb;
        Vec3 F, M;
        s->total(c.m, t0, F, M);
        // direction of the crest velocity of the 30 deg wave (positive u at the given time? use |F_h| and its angle)
        const double ang = std::atan2(F(1), F(0))*180.0/M_PI;
        const double Fh = std::hypot(F(0), F(1));
        const double zarm = std::hypot(M(0), M(1))/Fh;      // moment about the origin over the force: the zone's centre
        const bool dir = std::fabs(std::fabs(ang) - 30.0) < 1e-6 || std::fabs(std::fabs(ang) - 150.0) < 1e-6;
        check(std::fabs(Fh - Fpk) < 1e-9*Fpk && std::fabs(zarm - (1.0 - 0.5*lam)*eb) < 1e-9 && dir && std::fabs(F(2)) < 1e-6*Fh,
              fmt("peak %.4e N (rho R V^2 Cs(0) lambda eta_b %.4e), centre at %.4f m ((1 - lambda/2) eta_b), direction %.1f deg", Fh, Fpk, zarm, ang));
        // impulse
        double I = 0.0;
        const int n = 200000;
        const double T = slam_duration(cc)*R/V;
        for(int i=0; i<n; ++i)
        {
            s->total(c.m, t0 + (i + 0.5)*T/n, F, M);
            I += std::hypot(F(0), F(1))*T/n;
        }
        const double Iref = RHO*R*R*V*lam*eb*cs_integral(cc);
        s->total(c.m, t0 + T*1.0001, F, M);
        check(std::fabs(I - Iref) < 1e-4*Iref && F.norm()==0.0, fmt("impulse %.4e N s (rho R^2 V lambda eta_b int Cs %.4e), duration %.4f s", I, Iref, T));
    }
}

static void test_kinematic_trigger()
{
    std::cout<<"crest trigger, kinematic criterion, regular wave"<<std::endl;
    // deep water: u at the crest a w, c = g/w, B = a k
    const double H = 6.0, T = 10.0, a = 0.5*H, w = 2.0*M_PI/T, k = w*w/9.81;
    for(double Bc : {0.1, 0.2})
    {
        coral_case c;
        c.read_string(pile("waves regular H 6 T 10\n", fmt("model wo trigger kinematic B %.2f period 10 depth 1000 velocity crest", Bc)));
        member_slam *s = find_slam(c.m);
        s->initialize(c.m, 0.0);
        const double dt = 0.02;
        double tmax = 0.0, emax = 0.0, vmax = 0.0;
        int n0 = 0;
        for(int i=1; i<=(int)std::lround(105.0/dt); ++i)
        {
            const double t = i*dt;
            s->commit(c.m, t, dt);
            if(s->nevents() > n0)
            {
                n0 = s->nevents();
                const auto &e = s->events().back();
                // the crest at x = 0: cos(-w t) = 1, t = j T; the impact one step after it
                const double tc = std::round(e.t0/T)*T;
                tmax = std::max(tmax, std::fabs(e.t0 - tc - dt));
                emax = std::max(emax, std::fabs(e.eta_b - a));
                vmax = std::max(vmax, std::fabs(e.V - a*w));
            }
        }
        if(Bc < a*k)
            check(n0 == 10 && tmax < 1.5*dt && emax < 0.01*a && vmax < 0.01*a*w,
                  fmt("B = ak = %.3f >= %.2f: %.0f impacts in 10 periods, start <= %.3f s after crest + dt, eta_b and V within 1 %%", a*k, Bc, n0, tmax));
        else
            check(n0 == 0, fmt("B = ak = %.3f < %.2f: no impact", a*k, Bc));
    }
}

static void test_breaking_flag()
{
    std::cout<<"crest trigger, breaking flag of a sampled fluid"<<std::endl;
    coral_case c;
    c.read_string(pile("waves regular H 6 T 10\nwaves depth 30\nfluid_kinematics order 2\n", "model goda period 10"));
    member_slam *s = find_slam(c.m);
    s->initialize(c.m, 0.0);
    const double dt = 0.02;
    std::vector<int> imp;
    int n0 = 0;
    double V1 = 0.0;
    for(int i=0; i<=(int)std::lround(60.0/dt); ++i)
    {
        const double t = i*dt;
        std::vector<Vec3> pts;
        c.fluid_kin->collect(c.m, pts);
        std::vector<Vec3> u(pts.size());
        std::vector<double> eta(pts.size());
        std::vector<int> in(pts.size(), 1);
        for(size_t q=0; q<pts.size(); ++q)
        {
            eta[q] = c.waves->elevation(pts[q](0), pts[q](1), t);
            Vec3 p = pts[q], du;
            p(2) = std::min(p(2), eta[q]);
            double pr;
            c.waves->kinematics(p, t, u[q], du, pr);
            // breaking flagged on the rise of the crests at 20 s and 40 s
            const double tc = std::round(t/10.0)*10.0;
            if((tc == 20.0 || tc == 40.0) && t > tc - 1.0 && t <= tc)
                in[q] = 2;
        }
        c.fluid_kin->add_sample(t, u, eta, in);
        if(i == 0)
            continue;
        s->commit(c.m, t, dt);
        if(s->nevents() > n0) {n0 = s->nevents(); imp.push_back((int)std::lround(s->events().back().t0/10.0)); V1 = s->events().back().V;}
    }
    // the impact velocity: the celerity of T 10 s at the depth of the waves (30 m; the emulated samples have no bed)
    const double w = 2.0*M_PI/10.0, cel = w/wave_field::wavenumber(w, 30.0, 9.81);
    check(imp.size()==2 && imp[0]==2 && imp[1]==4,
          fmt("impacts only after the flagged crests: %.0f (at the crests of 20 s and 40 s)", (double)imp.size()));
    check(V1 > 0.0 && std::fabs(V1 - cel) < 1e-9*cel, fmt("impact velocity: the celerity %.4f m/s of T 10 s at the waves' depth (30 m)", V1));
}

static void test_state_keyword()
{
    std::cout<<"restart state, keyword"<<std::endl;
    coral_case c, d;
    const std::string in = pile("waves regular H 6 T 10\n", "model cw lambda 0.5 trigger time at 2.0 velocity 8 eta 5");
    c.read_string(in);
    d.read_string(in);
    member_slam *s = find_slam(c.m), *r = find_slam(d.m);
    s->initialize(c.m, 0.0);
    r->initialize(d.m, 0.0);
    for(int i=1; i<=100; ++i) s->commit(c.m, i*0.02, 0.02);
    std::vector<double> v;
    s->save_state(v);
    size_t k = 0;
    const bool ok = r->load_state(v, k) && k == v.size();
    Vec3 F1, F2, M;
    s->total(c.m, 2.1, F1, M);
    r->total(d.m, 2.1, F2, M);
    check(ok && F1.norm() > 0.0 && (F1 - F2).norm()==0.0 && r->nevents()==1, "save_state / load_state: the impact in progress continues");
    std::vector<std::string> n;
    std::vector<double> val;
    s->channel_names(n);
    s->channel_values(c.m, 2.1, val);
    coral_case b1, b2;
    b1.read_string(pile("waves regular H 6 T 10\n", "model wo period 12"));
    b2.read_string(pile("waves regular H 6 T 10\n", "model wo period 12 depth 23.5"));
    check(b1.sample_bed && !b2.sample_bed, "the celerity without depth asks the fluid for its bed");
    check(n.size()==9 && n[0]=="s_Fx" && n[8]=="s_n" && val.size()==9 && val[8]==1.0 && std::fabs(val[0] - F1(0)) < 1e-9*F1.norm(),
          "channels _Fx .. _Mz, _eta, _B, _n");
    auto err = [&](const std::string &slam, const char *msg)
    {
        try {coral_case e; e.read_string(pile("waves regular H 6 T 10\n", slam));}
        catch(const std::exception &x) {return std::string(x.what()).find(msg) != std::string::npos;}
        return false;
    };
    check(err("model fnv", "model must be") && err("model custom cs0 4", "'tau' is required") && err("trigger kinematic", "'period' is required")
          && err("model wo", "'period' is required") && err("trigger time", "'at <t1>") && err("lambda 1.5 period 10", "lambda"),
          "errors: unknown model, custom without tau, kinematic or celerity without period, time without at, lambda");
}

int main()
{
    test_coefficients();
    test_given_impact();
    test_kinematic_trigger();
    test_breaking_flag();
    test_state_keyword();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
