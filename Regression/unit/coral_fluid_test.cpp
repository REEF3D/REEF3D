// Standalone tests of the REEF3D::CORAL coupling to a resolved fluid (src/coral/hydro/coral_fluid_load,
// src/coral/core/coral_track): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         coral_fluid_test.cpp ../../src/coral/*/*.cpp (without reef3d/ and app/) -o coral_fluid_test
// Run:    ./coral_fluid_test
//
// Tests: body track (Hermite interpolation of position, velocity, rotation and angular velocity of an
// analytic motion between recorded states, exact at the recorded times); fluid load (extrapolation
// of the samples, added-mass force on the accelerations, channels); a heaving body driven by an
// emulated fluid solver that hands over the hydrostatic and damping loads of its own time levels plus
// the added mass (the coupling of coral/reef3d): period and decay against the damped oscillator
// (m + A) z'' + B z' + C z = 0: second order in the fluid time step; with quadratic extrapolation or the
// predictor (the load of the next fluid stage at a predicted pose as a provisional sample) the error of
// CORAL's own time integration; constant loads much worse (a lagged restoring force is a negative
// damping); the predictor stable at ten steps per period; a body 100 times lighter than its added
// mass stays stable; the fluid_body and fluid_coupling keywords.
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_track.h"
#include"coral/core/coral_rotation.h"
#include"coral/hydro/coral_fluid_load.h"
#include"coral/io/coral_input.h"
#include"coral/structure/coral_linear6.h"
#include<cmath>
#include<cstdio>
#include<iostream>
#include<memory>
#include<string>
using namespace coral;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}
static std::string fmt(const char *f, double a, double b=0.0, double c=0.0, double d=0.0)
{
    char s[256];
    std::snprintf(s, sizeof(s), f, a, b, c, d);
    return s;
}

// analytic motion: x(t), rotation about a fixed axis by th(t)
static const Vec3 axis = Vec3(1.0, 2.0, -0.5).normalized();
static Vec3 xa(double t) {return Vec3(std::sin(t), t*t, std::cos(2.0*t));}
static Vec3 va(double t) {return Vec3(std::cos(t), 2.0*t, -2.0*std::sin(2.0*t));}
static Vec3 aa(double t) {return Vec3(-std::sin(t), 2.0, -4.0*std::cos(2.0*t));}
static double th(double t) {return 0.3*std::sin(1.5*t);}
static double thd(double t) {return 0.45*std::cos(1.5*t);}
static double thdd(double t) {return -0.675*std::sin(1.5*t);}

static void set_state(model &m, double t)
{
    rigid_body &b = m.bodies[0];
    const Mat3 R0 = exp_so3(Vec3(0.2, -0.1, 0.4));
    b.x = xa(t); b.v = va(t); b.vd = aa(t);
    b.R = exp_so3(th(t)*axis)*R0;
    b.w = b.R.transpose()*(thd(t)*axis);
    b.wd = b.R.transpose()*(thdd(t)*axis);
}

static void test_track()
{
    std::cout<<"body track"<<std::endl;
    model m;
    rigid_body b;
    b.mass = 1.0; b.inertia = Mat3::Identity();
    m.add_body(b);
    body_track tr;
    tr.set_bodies({0});
    const double t0 = 1.0, h = 0.05;
    for(int k=0; k<=2; ++k)
    {
        set_state(m, t0 + k*h);
        tr.record(m, t0 + k*h);
    }
    double ex = 0.0, ev = 0.0, eR = 0.0, ew = 0.0, e0 = 0.0;
    for(double s : {0.1, 0.37, 0.5, 0.81, 1.0, 1.5, 1.93})
    {
        const double t = t0 + s*h;
        const body_snapshot q = tr.state(0, t);
        set_state(m, t);
        const rigid_body &r = m.bodies[0];
        ex = std::max(ex, (q.x - r.x).norm());
        ev = std::max(ev, (q.v - r.v).norm());
        eR = std::max(eR, log_so3(q.R*r.R.transpose()).norm());
        ew = std::max(ew, (q.w - r.R*r.w).norm());
        if(s == 1.0)
            e0 = (q.x - r.x).norm() + (q.v - r.v).norm() + log_so3(q.R*r.R.transpose()).norm();
    }
    // cubic Hermite: error h^4/384 max|f''''|
    check(ex < 2.0e-7 && ev < 1.0e-6, fmt("position and velocity between the records (max error %.1e m, %.1e m/s)", ex, ev));
    check(eR < 5.0e-8 && ew < 1.0e-7, fmt("rotation and angular velocity (max error %.1e rad, %.1e rad/s)", eR, ew));
    check(e0 < 1.0e-14, fmt("exact at a recorded time (%.1e)", e0));
    const body_snapshot qe = tr.state(0, t0 + 3.0*h);
    set_state(m, t0 + 2.0*h);
    check((qe.x - m.bodies[0].x).norm() < 1.0e-15, "clamped beyond the last record");
}

static void test_load()
{
    std::cout<<"fluid load"<<std::endl;
    model m;
    rigid_body b;
    b.name = "float"; b.mass = 2.0; b.inertia = Mat3::Identity();
    b.R = exp_so3(Vec3(0.0, 0.0, 0.5));
    m.add_body(b);
    fluid_load fl(m, {0});
    VecX F0(6), F1(6);
    F0 << 1, 2, 3, 4, 5, 6;
    F1 << 2, 2, 5, 4, 7, 6;
    MatX A = MatX::Zero(6, 6);
    A(2,2) = 10.0; A(4,4) = 3.0; A(2,4) = A(4,2) = 1.0;
    check(fl.force(0.3).norm() == 0.0, "no sample: no load");
    fl.set_sample(0.0, F0, A);
    check((fl.force(0.4) - F0).norm() == 0.0, "one sample: constant");
    fl.set_sample(0.1, F1);
    check((fl.force(0.15) - (F1 + 0.5*(F1 - F0))).norm() < 1.0e-12, "linear extrapolation between the fluid time levels");
    fl.order = 0;
    check((fl.force(0.15) - F1).norm() == 0.0, "constant with order 0");
    fl.order = 1;
    m.bodies[0].vd = Vec3(0.0, 0.0, 2.0);
    m.bodies[0].wd = m.bodies[0].R.transpose()*Vec3(0.0, 0.5, 0.0);
    VecX f = VecX::Zero(6);
    fl.eval(m, 0.1, f);
    const Vec3 Fz(2.0, 2.0, 5.0 - 10.0*2.0 - 1.0*0.5);
    const Vec3 Mg(4.0, 7.0 - 1.0*2.0 - 3.0*0.5, 6.0);
    check((f.head<3>() - Fz).norm() < 1.0e-12 && (f.tail<3>() - m.bodies[0].R.transpose()*Mg).norm() < 1.0e-12,
          "force and moment (body frame) with the added-mass part -A a");
    std::vector<std::string> names;
    fl.channel_names(names);
    check(names.size() == 6 && names[2] == "fluid_float_Fz", "channels fluid_<body>_Fx .. Mz");
}

// heaving body, emulated fluid: F = -C z - B z' at the fluid time levels, added mass A
static double heave_error(double mass, double A, double C, double B, double dt_fluid, int order, int sub, bool &stable,
                          bool implicit=false, bool predictor=false)
{
    model m;
    rigid_body b;
    b.name = "buoy"; b.mass = mass; b.inertia = Vec3(1, 1, 1).asDiagonal()*mass;
    m.add_body(b);
    m.gravity = Vec3::Zero();
    fluid_load &fl = static_cast<fluid_load&>(m.add_load(std::make_unique<fluid_load>(m, std::vector<int>{0})));
    fl.order = order;
    if(implicit)
    {
        Mat6 K = Mat6::Zero(), Cd = Mat6::Zero();
        K(2,2) = C; Cd(2,2) = B;
        m.add_load(std::make_unique<linear6>(m, 0, Vec3::Zero(), K, Cd));
    }
    const double z0 = 0.2;
    m.bodies[0].x(2) = z0;
    MatX Am = MatX::Zero(6, 6);
    Am(2,2) = A;
    auto sample = [&](double t)
    {
        VecX F = VecX::Zero(6);
        if(!implicit)
            F(2) = -C*m.bodies[0].x(2) - B*m.bodies[0].v(2);
        fl.set_sample(t, F, Am);
    };
    timeint_settings ts;
    ts.rho_inf = 1.0;
    timeint ti(ts);
    sample(0.0);
    ti.initialize(m, 0.0);
    const double M = mass + A, wn = std::sqrt(C/M), zeta = B/(2.0*std::sqrt(C*M)), wd = wn*std::sqrt(1.0 - zeta*zeta);
    const double T = 2.0*M_PI/wn;
    const int n = (int)std::lround(5.0*T/dt_fluid);
    double err = 0.0;
    stable = true;
    for(int i=0; i<n; ++i)
    {
        const double t = i*dt_fluid;
        if(i > 0) sample(t);
        if(predictor)
        {
            // the fluid's next stage with a predicted pose (Taylor), its load as a provisional sample
            const rigid_body &r = m.bodies[0];
            const double h = dt_fluid, zp = r.x(2) + h*r.v(2) + 0.5*h*h*r.vd(2), vp = r.v(2) + h*r.vd(2);
            VecX F = VecX::Zero(6);
            F(2) = -C*zp - B*vp;
            fl.set_sample(t + h, F, Am, true);
        }
        for(int k=0; k<sub; ++k)
            ti.step(m, t + k*dt_fluid/sub, dt_fluid/sub);
        const double t1 = t + dt_fluid;
        const double ze = z0*std::exp(-zeta*wn*t1)*(std::cos(wd*t1) + zeta*wn/wd*std::sin(wd*t1));
        err = std::max(err, std::fabs(m.bodies[0].x(2) - ze));
        if(!(std::fabs(m.bodies[0].x(2)) < 10.0*z0)) {stable = false; break;}
    }
    return err/z0;
}

static void test_coupled_heave()
{
    std::cout<<"heaving body driven by an emulated fluid solver"<<std::endl;
    const double mass = 1000.0, A = 600.0, C = 4.0e4, B = 800.0;
    const double T = 2.0*M_PI*std::sqrt((mass + A)/C);
    bool st;
    const double ei = heave_error(mass, A, C, B, T/100.0, 1, 1, st, true);
    const double e1a = heave_error(mass, A, C, B, T/50.0, 1, 1, st);
    const double e1b = heave_error(mass, A, C, B, T/100.0, 1, 1, st);
    const double e2b = heave_error(mass, A, C, B, T/100.0, 2, 1, st);
    const double epa = heave_error(mass, A, C, B, T/50.0, 2, 1, st, false, true);
    const double epb = heave_error(mass, A, C, B, T/100.0, 2, 1, st, false, true);
    const double e0b = heave_error(mass, A, C, B, T/100.0, 0, 1, st);
    // the time integration alone (spring and damper in CORAL) has the period error (w h)^2/12 of the
    // trapezoidal rule (rho_inf = 1); linear extrapolation adds about four times that
    check(e1b < 6.0*ei && e1a/e1b > 3.5 && e1a/e1b < 4.5,
          fmt("linear extrapolation: max error %.2e of z0 over 5 periods at T/100 (CORAL alone %.2e), ratio %.2f from T/50", e1b, ei, e1a/e1b));
    check(e2b < 1.5*ei, fmt("quadratic extrapolation: max error %.2e of z0", e2b));
    check(epb < 1.5*ei && epa/epb > 3.5 && epa/epb < 4.5, fmt("predictor (provisional load of the next stage): max error %.2e of z0, ratio %.2f from T/50", epb, epa/epb));
    check(e0b > 10.0*e1b, fmt("constant loads lag the restoring force by a step (negative damping): max error %.2e of z0", e0b));
    const double es = heave_error(mass, A, C, B, T/100.0, 2, 4, st, false, true);
    check(es < epb, fmt("four CORAL sub-steps per fluid step: max error %.2e of z0", es));
    bool s1, s2, s3;
    const double e10 = heave_error(mass, A, C, B, T/10.0, 2, 1, s1, false, true);
    heave_error(6.0, A, C, B, T/10.0, 2, 1, s2, false, true);
    heave_error(6.0, A, C, B, T/10.0, 1, 1, s3);
    check(s1 && s2 && !s3, fmt("ten fluid steps per period: stable with the predictor (max error %.2f of z0, also the light body), the light body unstable with linear extrapolation", e10));
    // light body: m = A/100 (the added mass explicit in the loads would be unstable)
    const double el = heave_error(6.0, A, C, B, T/100.0, 2, 1, st, false, true);
    check(st && el < 3.0*epb, fmt("body 100 times lighter than its added mass: stable, max error %.2e of z0", el));
}

static void test_input()
{
    std::cout<<"fluid_body and fluid_coupling keywords"<<std::endl;
    coral_case c;
    c.read_string("body float mass 100 cog 0 0 -0.5\nbody spar mass 200 cog 0 0 -5\n"
                  "fluid_body spar n 1\nfluid_body float n 0 own_mass\nfluid_coupling order 1 no_predictor\n");
    check(c.fluid != nullptr && c.fluid->nbody() == 2 && c.fluid->bodies()[0] == 1 && c.fluid->order == 1 && !c.fluid_predictor
          && c.fluid_bodies[1].n == 0 && c.fluid_bodies[1].own_mass && !c.fluid_bodies[0].own_mass,
          "two fluid bodies, REEF3D numbers, own_mass, order, no_predictor");
    bool err = false;
    try {coral_case d; d.read_string("body a mass 1\nfluid_body a\nfluid_body a\n");}
    catch(const std::exception &e) {err = std::string(e.what()).find("twice") != std::string::npos;}
    check(err, "a body mapped twice is an error");
}

int main()
{
    test_track();
    test_load();
    test_coupled_heave();
    test_input();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
