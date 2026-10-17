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
// mass stays stable; the fluid_body and fluid_coupling keywords. A fluid without the instantaneous added
// mass (NHFLOW: the reaction of the body acceleration it saw is in its loads): the artificial added mass
// of CORAL stabilises a body ten times lighter than its added mass; As = A is as accurate as the implicit
// added mass, other values are stable with an error growing with |A - As|.
// Sampled kinematics (src/coral/hydro/coral_sampled_kinematics): the points of Morison members on frame
// members and bodies and of a dynamic mooring line (kinematics_points); an emulated fluid (linear waves)
// sampled there once per fluid step: velocity, acceleration (time derivative of the polynomial) and
// elevation against the waves themselves, the Morison force of a monopile (frame member, fluid_kinematics
// in the input) against the same model with CORAL's waves (second order in the fluid step); points outside
// the fluid domain and requests before the first sample take the background waves.
// One-way loads: fluid_body ... one_way (load at the REEF3D CoG on another body: force and moment),
// fluid_member (strip forces of a resolved hull along a frame member: element and weight of each strip,
// total force and moment, the static deflection of a cantilever under a uniform strip load against
// q L^4/(8 EI), excitation without the first sample, time polynomial).
// Record and replay (src/coral/io/coral_coupling): a coupling driven by an emulated fluid (kinematics,
// one-way body loads, strip loads, steps) writes its journal; the replay on a fresh model repeats the
// time series bit for bit.
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_track.h"
#include"coral/core/coral_rotation.h"
#include"coral/hydro/coral_fluid_load.h"
#include"coral/io/coral_input.h"
#include"coral/structure/coral_linear6.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include"coral/hydro/coral_member_fluid_load.h"
#include"coral/io/coral_coupling.h"
#include<filesystem>
#include<fstream>
#include<sstream>
#include"coral/core/coral_static.h"
#include"coral/structure/coral_beam.h"
#include"coral/hydro/coral_member_morison.h"
#include"coral/hydro/coral_morison.h"
#include"coral/mooring/coral_mooring.h"
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

// heaving body, emulated fluid without the instantaneous added mass (NHFLOW): its loads contain the
// reaction -A a_f of the body acceleration a_f it saw (the mean of its step), F = -C z - B z' - A a_f;
// CORAL's artificial added
// mass As makes the staggered coupling implicit with compensation (fluid_load::set_artificial_added_mass)
static double heave_error_lagged(double mass, double A, double As, double C, double B, double dt_fluid, bool &stable,
                                 bool predictor=true, int order=2)
{
    model m;
    rigid_body b;
    b.name = "buoy"; b.mass = mass; b.inertia = Vec3(1, 1, 1).asDiagonal()*mass;
    m.add_body(b);
    m.gravity = Vec3::Zero();
    fluid_load &fl = static_cast<fluid_load&>(m.add_load(std::make_unique<fluid_load>(m, std::vector<int>{0})));
    fl.order = order;
    Mat6 Ma = Mat6::Zero();
    Ma(2,2) = As;
    if(As > 0.0)
        fl.set_artificial_added_mass(0, Ma);
    const double z0 = 0.2;
    m.bodies[0].x(2) = z0;
    // the accelerations of the body when the fluid sampled (as the coupling passes them)
    auto sample = [&](double t, double z, double v, double a, bool prov)
    {
        VecX F = VecX::Zero(6), acc = VecX::Zero(6);
        F(2) = -C*z - B*v - A*a;
        acc(2) = a;
        fl.set_sample(t, F, MatX(), prov, acc);
    };
    timeint_settings ts;
    ts.rho_inf = 1.0;
    timeint ti(ts);
    // start at rest: the consistent acceleration of the free decay
    const double a0 = -C*z0/(mass + A);
    m.bodies[0].vd(2) = a0;
    sample(0.0, z0, 0.0, a0, false);
    ti.initialize(m, 0.0);
    const double M = mass + A, wn = std::sqrt(C/M), zeta = B/(2.0*std::sqrt(C*M)), wd = wn*std::sqrt(1.0 - zeta*zeta);
    const double T = 2.0*M_PI/wn;
    const int n = (int)std::lround(5.0*T/dt_fluid);
    double err = 0.0;
    stable = true;
    for(int i=0; i<n; ++i)
    {
        const double t = i*dt_fluid, h = dt_fluid;
        const rigid_body &r = m.bodies[0];
        // predictor: the fluid's stage at the predicted pose (Taylor, constant acceleration)
        if(predictor)
            sample(t + h, r.x(2) + h*r.v(2) + 0.5*h*h*r.vd(2), r.v(2) + h*r.vd(2), r.vd(2), true);
        const double v0 = r.v(2);
        ti.step(m, t, h);
        // the fluid at the end of the step, at CORAL's motion: it reacts to the mean acceleration of its
        // step (the change of the imposed body velocity, as NHFLOW's direct forcing; the coupling passes
        // the same mean)
        sample(t + h, r.x(2), r.v(2), (r.v(2) - v0)/h, false);
        const double t1 = t + h;
        const double ze = z0*std::exp(-zeta*wn*t1)*(std::cos(wd*t1) + zeta*wn/wd*std::sin(wd*t1));
        err = std::max(err, std::fabs(r.x(2) - ze));
        if(!(std::fabs(r.x(2)) < 10.0*z0)) {stable = false; break;}
    }
    return err/z0;
}

static void test_artificial_added_mass()
{
    std::cout<<"fluid without its own added mass: artificial added mass"<<std::endl;
    const double A = 600.0, C = 4.0e4, B = 800.0;
    const double T = 2.0*M_PI*std::sqrt((1000.0 + A)/C);
    bool st, s0, s1, s2, s3;
    // reference: the same body with the added mass implicit (as FNPF hands it over)
    const double eref = heave_error(1000.0, A, C, B, T/100.0, 2, 1, st, false, true);
    const double eh = heave_error_lagged(1000.0, A, A, C, B, T/100.0, s1);
    const double eh0 = heave_error_lagged(1000.0, A, 0.0, C, B, T/100.0, s0);
    check(s1 && eh < 3.0*eref, fmt("m = 1000, A = 600: As = A max error %.2e of z0 (added mass implicit %.2e; As = 0: %.2e)", eh, eref, s0 ? eh0 : 1.0/0.0));
    // light body m = A/10: explicit added mass unstable, As stabilises, the result does not depend on As
    const double T2 = 2.0*M_PI*std::sqrt((60.0 + A)/C);
    heave_error_lagged(60.0, A, 0.0, C, B, T2/100.0, s0);
    const double ea = heave_error_lagged(60.0, A, A, C, B, T2/100.0, s1);
    const double eb = heave_error_lagged(60.0, A, 0.6*A, C, B, T2/100.0, s2);
    const double ec = heave_error_lagged(60.0, A, 2.0*A, C, B, T2/100.0, s3);
    check(!s0, "light body (m = A/10) without artificial added mass: unstable");
    // the provisional sample carries -(A - As) a^n explicitly: As = A is as accurate as the implicit
    // added mass, other values stay stable with an error that grows with |A - As|
    check(s1 && s2 && s3 && ea < 1.5*heave_error(60.0, A, C, B, T2/100.0, 2, 1, st, false, true) && eb > ea && ec > ea && eb < 0.25 && ec < 0.25,
          fmt("light body with As = A, 0.6 A, 2 A: stable, max error %.2e %.2e %.2e of z0 over 5 periods", ea, eb, ec));
    const double ea2 = heave_error_lagged(60.0, A, A, C, B, T2/200.0, st);
    check(st && ea/ea2 > 3.0, fmt("light body, As = A: second order in the fluid step (error ratio %.2f from T/100 to T/200)", ea/ea2));
    const double en = heave_error_lagged(60.0, A, A, C, B, T2/100.0, st, false, 2);
    check(st && en < 2.0*ea, fmt("light body, As = A, without the predictor (quadratic extrapolation): stable, max error %.2e of z0", en));
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

// one fluid step: sample the emulated fluid (linear waves w) at the points of the sampled source;
// inside: points with x < xmax only
static void emulate_sample(const model &m, sampled_kinematics &sk, const wave_field &w, double t, double xmax=1.0e30)
{
    std::vector<Vec3> pts;
    sk.collect(m, pts);
    std::vector<Vec3> u(pts.size());
    std::vector<double> eta(pts.size());
    std::vector<int> in(pts.size());
    for(size_t q=0; q<pts.size(); ++q)
    {
        Vec3 du; double pr;
        w.kinematics(pts[q], t, u[q], du, pr);
        eta[q] = w.elevation(pts[q](0), pts[q](1), t);
        in[q] = pts[q](0) < xmax ? 1 : 0;
    }
    sk.add_sample(t, u, eta, in);
}

static member_morison *find_member_morison(const model &m)
{
    for(const auto &l : m.load_list)
        if(auto *p = dynamic_cast<member_morison*>(l.get()))
            return p;
    return nullptr;
}

static void test_sampled()
{
    std::cout<<"sampled kinematics of a resolved fluid"<<std::endl;
    // monopile D 6 m in 20 m of water, regular wave H 4 m, T 8 s; the same model once with CORAL's
    // waves and once with fluid_kinematics (the waves of coral.txt are then the background)
    const std::string base =
        "rho 1025\nwaves regular H 4 T 8\nwaves depth 20\n"
        "section pile tube D 6 t 0.06 E 2.1e11 G 8.08e10 rho 7850\n"
        "member mp p1 0 0 -20 p2 0 0 10 n 10 section pile orient 1 0 0 end1 ground\n"
        "morison wave member mp Cd 1.0 Ca 1.0 n 4\n";
    coral_case ca, cb;
    ca.read_string(base);
    cb.read_string(base + "fluid_kinematics order 2\n");
    member_morison *ma = find_member_morison(ca.m), *mb = find_member_morison(cb.m);
    check(cb.fluid_kin != nullptr && ma != nullptr && mb != nullptr, "fluid_kinematics: the Morison member takes the sampled source");
    std::vector<Vec3> pts;
    cb.fluid_kin->collect(cb.m, pts);
    std::vector<Vec3> pa;
    mb->kinematics_points(cb.m, ca.waves.get(), pa);
    check(pts.size() == 40 && pa.empty() && std::fabs(pts[0](2) - (-20.0 + 3.0/8.0)) < 1.0e-12,
          fmt("points: %.0f strip centres of the member (and none for another source)", (double)pts.size()));
    // before the first sample: the background waves
    Vec3 Fa, Ma, Fb, Mb;
    ma->total(ca.m, 1.3, Fa, Ma);
    mb->total(cb.m, 1.3, Fb, Mb);
    check((Fa - Fb).norm() < 1.0e-9*Fa.norm(), "before the first sample: the background waves");
    // fluid steps of T/40 and T/80: after each sample the force over the last interval
    const wave_field &w = *ca.waves;
    double err[2];
    for(int r=0; r<2; ++r)
    {
        const double dt = 8.0/(40.0*(r + 1));
        sampled_kinematics &sk = *cb.fluid_kin;
        double e = 0.0, fmax = 0.0;
        for(int k=0; k<=(int)std::lround(16.0/dt); ++k)
        {
            const double tk = k*dt;
            emulate_sample(cb.m, sk, w, tk);
            if(k < 2)
                continue;
            for(double s : {0.25, 0.5, 0.75, 1.0})
            {
                const double t = tk - (1.0 - s)*dt;
                ma->total(ca.m, t, Fa, Ma);
                mb->total(cb.m, t, Fb, Mb);
                e = std::max(e, (Fa - Fb).norm());
                fmax = std::max(fmax, Fa.norm());
            }
        }
        err[r] = e/fmax;
    }
    check(err[1] < 1.0e-3 && err[0]/err[1] > 3.5, fmt("monopile Morison force: max error %.1e (T/80) of the force amplitude, ratio %.1f from T/40 (second order: the acceleration is the time derivative)", err[1], err[0]/err[1]));
    // velocity, acceleration and elevation at one point between the samples
    {
        Vec3 u, du, uw, duw; double pr;
        const double t = 16.0 - 0.37*0.1;
        const Vec3 p = pts[30];
        mb->kinematics_points(cb.m, cb.fluid_kin.get(), pa);
        cb.fluid_kin->kinematics_at(mb, 30, p, t, u, du, pr);
        w.kinematics(p, t, uw, duw, pr);
        const double eu = (u - uw).norm()/2.0, ed = (du - duw).norm()/2.0;
        const double ee = std::fabs(cb.fluid_kin->elevation_at(mb, 30, p, t) - w.elevation(p(0), p(1), t));
        check(eu < 1.0e-4 && ed < 2.0e-3 && ee < 1.0e-4, fmt("velocity, acceleration, elevation between the samples (errors %.1e, %.1e, %.1e)", eu, ed, ee));
    }
    // points outside the fluid domain: background; without background: still water
    {
        model m;
        rigid_body b;
        b.name = "buoy"; b.mass = 100.0; b.inertia = Mat3::Identity();
        m.add_body(b);
        auto bg = std::make_shared<wave_field>(20.0);
        bg->add_regular(2.0, 6.0);
        auto sk = std::make_shared<sampled_kinematics>(bg);
        morison_member mm;
        mm.p1 = Vec3(-5.0, 0.0, -6.0); mm.p2 = Vec3(5.0, 0.0, -6.0); mm.D = 0.5; mm.n = 10;
        m.add_load(std::make_unique<morison_load>(m, 0, std::vector<morison_member>{mm}, sk, 1025.0));
        mooring_end ea, eb;
        ea.body = -1; ea.p = Vec3(-30.0, 0.0, -20.0);
        eb.body = 0; eb.p = Vec3(0.0, 0.0, 0.0);
        mooring_props mp;
        mp.L = 40.0; mp.d = 0.1; mp.mass = 20.0; mp.EA = 1.0e8;
        m.add_load(std::make_unique<mooring_line>(m, ea, eb, mp, 8, sk));
        std::vector<Vec3> q;
        sk->collect(m, q);
        check(q.size() == 10 + 9, fmt("points of a Morison member on a body and of a dynamic line: %.0f", (double)q.size()));
        // emulated fluid = other waves, in the domain only for x < 0
        wave_field wf(20.0);
        wf.add_regular(3.0, 9.0);
        for(int k=0; k<3; ++k)
            emulate_sample(m, *sk, wf, k*0.1, 0.0);
        const void *owner = m.load_list[0].get();
        Vec3 u, du, ue, due, ub, dub; double pr;
        sk->kinematics_at(owner, 0, q[0], 0.15, u, du, pr);      // x = -4.5: inside
        wf.kinematics(q[0], 0.15, ue, due, pr);
        sk->kinematics_at(owner, 9, q[9], 0.15, ub, dub, pr);    // x = 4.5: outside
        Vec3 ubg, dubg;
        bg->kinematics(q[9], 0.15, ubg, dubg, pr);
        check((u - ue).norm() < 1.0e-4 && (ub - ubg).norm() < 1.0e-12, "inside the fluid domain the samples, outside the background waves");
        sampled_kinematics none;
        none.kinematics_at(owner, 0, q[0], 0.1, u, du, pr);
        check(u.norm() == 0.0 && none.elevation(0.0, 0.0, 0.1) == 0.0, "without samples and background: still water");
    }
}

static void test_one_way()
{
    std::cout<<"one-way loads: fluid_body ... one_way, fluid_member"<<std::endl;
    // a rigid body loaded at a point that is not its CoG (one-way load of a fluid body)
    {
        model m;
        rigid_body b;
        b.name = "tp"; b.mass = 10.0; b.inertia = Mat3::Identity(); b.x = Vec3(1.0, 0.0, 0.0);
        m.add_body(b);
        fluid_load fl(m, {0});
        fl.set_attachment(m, 0, Vec3(1.0, 0.0, 5.0));
        VecX F(6);
        F << 100.0, 0.0, 0.0, 0.0, 7.0, 0.0;
        fl.set_sample(0.0, F);
        VecX f = VecX::Zero(6);
        fl.eval(m, 0.0, f);
        // force 100 N at 5 m above the CoG: moment about the CoG 500 N m (y) plus the 7 N m of the load
        check(std::fabs(f(0) - 100.0) < 1.0e-12 && std::fabs(f(4) - 507.0) < 1.0e-9, fmt("load at the REEF3D CoG: My about the CORAL CoG %.1f N m (507)", f(4)));
    }
    // cantilever of frame elements (fixed at z = -20), a uniform horizontal load of a resolved hull in strips
    const std::string base =
        "gravity off\nsection pile tube D 6 t 0.06 E 2.1e11 G 8.08e10 rho 7850\n"
        "member mp p1 0 0 -20 p2 0 0 10 n 30 section pile orient 1 0 0 end1 ground\n";
    coral_case c;
    c.read_string(base + "fluid_member hull member mp n 0 strips 60\n");
    check(c.fluid_members.size() == 1 && c.fluid_members[0].n == 0 && c.fluid_members[0].load->nstrips() == 60
          && std::fabs(c.fluid_members[0].load->length() - 30.0) < 1.0e-12, "fluid_member: 60 strips along the 30 m member, 6DOF body 0");
    member_fluid_load *l = c.fluid_members[0].load;
    const double q = 2.0e4;                     // N/m
    VecX F = VecX::Zero(3*60);
    for(int k=0; k<60; ++k) F(3*k) = q*0.5;     // strips of 0.5 m
    l->set_sample(0.0, F);
    Vec3 Ft, Mt;
    l->total(c.m, 0.0, Ft, Mt);
    // moment about the origin: the load acts between z = -20 and 10, centre -5
    check(std::fabs(Ft(0) - 30.0*q) < 1.0e-6 && std::fabs(Mt(1) - (-5.0)*30.0*q) < 1.0e-3,
          fmt("total %.0f N, moment %.0f N m about the origin", Ft(0), Mt(1)));
    const static_result sr = static_solve(c.m, 0.0);
    const int top = c.m.find_body("mp:30");
    const beam_section s = beam_section::tube(6.0, 0.06, 2.1e11, 8.08e10, 7850.0);
    const double L = 30.0, wEB = q*L*L*L*L/(8.0*s.EIy), wS = q*L*L/(2.0*s.GAy);
    const double w = c.m.body(top).x(0);
    check(sr.converged && std::fabs(w - (wEB + wS)) < 2.0e-3*(wEB + wS), fmt("cantilever tip under the strip load %.5f m vs q L^4/8EI + q L^2/2GA = %.5f m", w, wEB + wS));
    // excitation and time polynomial
    coral_case c2;
    c2.read_string(base + "fluid_member hull member mp strips 2 excitation order 1\n");
    member_fluid_load *l2 = c2.fluid_members[0].load;
    VecX a(6), b(6);
    a << 1, 2, 3, 4, 5, 6;
    b << 2, 2, 5, 4, 7, 9;
    l2->set_sample(0.0, a);
    l2->set_sample(1.0, b);
    const VecX g = l2->forces(1.5);
    check(l2->excitation && l2->order == 1 && std::fabs(g(0) - 1.5) < 1.0e-12 && std::fabs(g(5) - 4.5) < 1.0e-12,
          "excitation: the first sample removed; linear in time");
}

static void test_replay()
{
    std::cout<<"record and replay of the fluid input"<<std::endl;
    const std::string text =
        "rho 1025\ntime dt 0.01 end 10\noutput dt 0.05 file cpl.dat\n"
        "waves regular H 2 T 6\nwaves depth 20\nfluid_kinematics\nfluid_coupling record cpl.jrn origin 0 0 20\n"
        "section pile tube D 4 t 0.04 E 2.1e11 G 8.08e10 rho 7850\n"
        "member mp p1 0 0 -20 p2 0 0 10 n 10 section pile orient 1 0 0 end1 ground\n"
        "morison wave member mp Cd 1 Ca 1 n 2\n"
        "fluid_member hull member mp n 0 strips 20\n"
        "body buoy mass 1000 inertia 100 100 100 cog 20 0 0\n"
        "joint heave prismatic ground buoy point 20 0 0 axis 0 0 1\n"
        "linear6 spring body buoy point 20 0 0 k 0 0 3e4 0 0 0\n"
        "fluid_body buoy n 1 one_way\n";
    namespace fs = std::filesystem;
    const fs::path d1 = "coral_fluid_test_rec", d2 = "coral_fluid_test_replay";
    fs::remove_all(d1); fs::remove_all(d2);
    fs::create_directories(d1); fs::create_directories(d2);
    std::ostringstream log;
    {
        coral_case c;
        c.read_string(text);
        coupling cp(c, log, true);
        cp.record(d1.string());
        wave_field w(20.0);
        w.add_regular(2.0, 6.0);
        auto sample = [&](double t)
        {
            std::vector<Vec3> pts;
            cp.kinematics_points(pts);
            std::vector<Vec3> u(pts.size());
            std::vector<double> eta(pts.size());
            std::vector<int> in(pts.size(), 1);
            for(size_t q=0; q<pts.size(); ++q)
            {
                Vec3 du; double pr;
                const Vec3 pc = pts[q] - Vec3(0.0, 0.0, 20.0);
                w.kinematics(pc, t, u[q], du, pr);
                eta[q] = 20.0 + w.elevation(pc(0), pc(1), t);
            }
            cp.kinematics(t, u, eta, in);
        };
        sample(0.0);
        cp.ini(0.0, d1.string());
        Vec3 a, e; double L; int ns;
        cp.connect(0, 5.0e5, Vec3(0.0, 0.0, 5.0), Mat3::Identity(), Mat3::Identity());
        cp.connect(1, 1.0e3, Vec3(20.0, 0.0, 20.0), Mat3::Identity(), Mat3::Identity());
        const bool st = cp.strips(0, a, e, L, ns);
        check(st && ns == 20 && std::fabs(L - 30.0) < 1.0e-12 && std::fabs(a(2)) < 1.0e-12 && cp.one_way(1) && cp.one_way(0),
              "strips of the fluid member in the fluid frame (origin 0 0 20); one-way bodies");
        const double dt = 0.07;
        for(int k=0; k<40; ++k)
        {
            const double t = k*dt;
            // FNPF-like: loads of the start of the step, then the fluid step, then the CORAL step
            VecX Fs = VecX::Zero(3*20);
            for(int i=0; i<20; ++i) Fs(3*i) = 1.0e4*std::sin(1.3*t + 0.1*i);
            cp.fluid_strips(0, t, Fs);
            std::vector<Vec6> F(2, Vec6::Zero());
            F[1](2) = 5.0e3*std::sin(1.1*t);
            cp.fluid_loads(t, F, {}, false);
            sample(t + dt);
            cp.step_to(t + dt);
        }
    }
    {
        coral_case c;
        c.read_string(text);
        std::ostringstream log2;
        const long nev = coupling::replay(c, (d1/"cpl.jrn").string(), d2.string(), log2);
        std::ifstream a1(d1/"cpl.dat"), a2(d2/"cpl.dat");
        std::stringstream s1, s2;
        s1<<a1.rdbuf(); s2<<a2.rdbuf();
        const bool same = s1.str() == s2.str() && s1.str().size() > 1000;
        check(nev > 100 && same, fmt("replay of %.0f events: time series identical (%.0f bytes)", (double)nev, (double)s1.str().size()));
    }
    fs::remove_all(d1); fs::remove_all(d2);
}

int main()
{
    test_track();
    test_load();
    test_coupled_heave();
    test_artificial_added_mass();
    test_input();
    test_sampled();
    test_one_way();
    test_replay();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
