// Standalone tests of the REEF3D::CORAL mooring lines (src/coral/mooring): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         coral_mooring_test.cpp ../../src/coral/core/*.cpp ../../src/coral/structure/*.cpp
//         ../../src/coral/hydro/*.cpp ../../src/coral/mooring/*.cpp -o coral_mooring_test
// Run:    ./coral_mooring_test
//
// Tests: elastic catenary (fairlead forces from the end position and back, hanging free and with a
// seabed portion and friction, slack line with length to spare, profile ending at the fairlead,
// anchor forces), quasi-static line on a
// body (force direction, symmetric spread), lumped-mass line relaxed between fixed ends against the
// catenary (fairlead tension, seabed length), surge decay of a moored body with the dynamic line
// against the quasi-static line (period), drag of a taut line in a current, a line between two bodies
// (no seabed by default: quasi-static and dynamic line hang as the free catenary) and an anchor above the
// seabed (hanging free), wave kinematics once per step against every Runge-Kutta stage, the end nodes'
// inertia on the bodies (rigid motion of a line) and the step independence of the fairlead tension
// (Hermite end motion), buoyancy of a line through the free surface, seabed friction along and across,
// the stop of an unstable line, the start from the static equilibrium of the segments, mooring systems with
// free points (clump weight, buoy, bridle: static balance; pendulum period; clump on the seabed), line
// failure at a time and at a tension, restart after a failure, seabed on a slope, a grid and the bed of a
// sampled fluid, bending stiffness (cantilever: deflection, clamping moment; frequency of a taut line).
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_gravity.h"
#include"coral/structure/coral_linear6.h"
#include"coral/mooring/coral_catenary.h"
#include"coral/mooring/coral_mooring.h"
#include"coral/hydro/coral_waves.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include<iostream>
#include<stdexcept>
#include<ctime>
#include<string>
#include<cstdio>
#include<cmath>
#include<memory>

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

static void test_catenary()
{
    std::cout<<"elastic catenary"<<std::endl;
    catenary_line c;
    c.L = 100.0; c.w = 500.0; c.EA = 1.0e8;
    // hanging free (no seabed)
    c.CB = -1.0;
    {
        const double H0 = 2.0e4, V0 = 6.0e4;
        double x, z;
        catenary_ends(c, H0, V0, x, z);
        double H = -1.0, V = 0.0;
        const bool ok = catenary_solve(c, x, z, H, V);
        check(ok && std::fabs(H - H0) < 1.0e-6*H0 && std::fabs(V - V0) < 1.0e-6*V0, fmt("hanging free: HF %.3f kN, VF %.3f kN back from the fairlead position", H/1e3, V/1e3));
        double xp, zp;
        catenary_profile(c, H, V, c.L, xp, zp);
        check(std::fabs(xp - x) < 1.0e-8 && std::fabs(zp - z) < 1.0e-8, "hanging free: profile ends at the fairlead");
        double HA, VA;
        catenary_anchor(c, H, V, HA, VA);
        check(std::fabs(HA - H) < 1.0e-9 && std::fabs(VA - (V - c.w*c.L)) < 1.0e-6, "hanging free: anchor HA = HF, VA = VF - w L");
    }
    // seabed portion with friction
    c.CB = 0.5;
    {
        const double H0 = 1.0e4, V0 = 3.0e4;    // 40 m on the seabed
        double x, z;
        catenary_ends(c, H0, V0, x, z);
        double H = -1.0, V = 0.0;
        const bool ok = catenary_solve(c, x, z, H, V);
        check(ok && std::fabs(H - H0) < 1.0e-6*H0 && std::fabs(V - V0) < 1.0e-6*V0, fmt("seabed and friction: HF %.3f kN, VF %.3f kN back", H/1e3, V/1e3));
        check(std::fabs(catenary_seabed_length(c, H, V) - 40.0) < 1.0e-6, fmt("seabed length %.4f m", catenary_seabed_length(c, H, V)));
        double xp, zp;
        catenary_profile(c, H, V, c.L, xp, zp);
        check(std::fabs(xp - x) < 1.0e-8 && std::fabs(zp - z) < 1.0e-8, fmt("seabed: profile ends at the fairlead (%.1e)", std::fabs(xp - x) + std::fabs(zp - z)));
        double HA, VA;
        catenary_anchor(c, H, V, HA, VA);
        check(std::fabs(HA - std::max(0.0, H0 - 0.5*500.0*40.0)) < 1.0e-6 && VA==0.0, fmt("anchor: HA = HF - CB w LB = %.1f N", HA));
    }
    // slack line: length to spare (xF + zF < L): no horizontal force, hanging vertically
    {
        catenary_line s = c;
        s.CB = 0.0;
        double H = -1.0, V = 0.0;
        const bool ok = catenary_solve(s, 80.0, 15.0, H, V);
        check(ok && H==0.0 && std::fabs(V - s.w*15.0) < 1.0e-3*s.w*15.0, fmt("slack line: HF = 0, VF = w zF (%.1f N vs %.1f N)", V, s.w*15.0));
        double H2 = -1.0, V2 = 0.0;
        const bool ok2 = catenary_solve(s, 85.5, 15.0, H2, V2);
        check(ok2 && H2 > 0.0 && H2 < 0.05*s.w*s.L, fmt("just taut (xF + zF = 1.005 L): small HF %.2f N", H2));
    }
    // inextensible limit: the classic catenary z = H/w (cosh(w x/H) - 1) from the touchdown point
    c.EA = 1.0e14;
    c.CB = 0.0;
    {
        const double H0 = 1.0e4, V0 = 3.0e4;
        double x, z;
        catenary_ends(c, H0, V0, x, z);
        const double LB = c.L - V0/c.w;
        const double xs = x - LB;
        check(std::fabs(z - H0/c.w*(std::cosh(c.w*xs/H0) - 1.0)) < 1.0e-6, "inextensible: z = H/w (cosh(w x/H) - 1)");
    }
}

static mooring_props chain()
{
    mooring_props p;
    p.L = 120.0;
    p.d = 0.1;
    p.mass = 50.0;
    p.EA = 1.0e8;
    p.BA = 2.0e5;
    p.CB = 0.0;
    p.Cdn = 1.2;
    p.Cdt = 0.05;
    p.Can = 1.0;
    return p;
}

static void test_quasi_static()
{
    std::cout<<"quasi-static line on a body"<<std::endl;
    model m;
    rigid_body b;
    b.mass = 1.0e5;
    b.inertia = Vec3(1e7, 1e7, 1e7).asDiagonal();
    m.add_body(b);
    const mooring_props p = chain();
    VecX f = VecX::Zero(6);
    for(int k=0; k<3; ++k)
    {
        const double a = 2.0*M_PI*k/3.0;
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(100.0*std::cos(a), 100.0*std::sin(a), -50.0);
        fl.body = 0;  fl.p = Vec3(5.0*std::cos(a), 5.0*std::sin(a), -5.0);
        mooring_qs l(m, an, fl, p);
        l.eval(m, 0.0, f);
    }
    check(std::fabs(f(0)) < 1.0e-6 && std::fabs(f(1)) < 1.0e-6 && f(2) < 0.0, fmt("three symmetric lines: no horizontal force, downward %.1f kN", -f(2)/1e3));
    // offset: restoring
    m.bodies[0].x = Vec3(3.0, 0.0, 0.0);
    VecX f2 = VecX::Zero(6);
    for(int k=0; k<3; ++k)
    {
        const double a = 2.0*M_PI*k/3.0;
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(100.0*std::cos(a), 100.0*std::sin(a), -50.0);
        fl.body = 0;  fl.p = Vec3(5.0*std::cos(a), 5.0*std::sin(a), -5.0);
        model m0;
        m0.add_body(b);
        mooring_qs l(m0, an, fl, p);
        l.eval(m, 0.0, f2);
    }
    check(f2(0) < 0.0 && std::fabs(f2(1)) < 1.0e-6, fmt("surge offset 3 m: restoring force %.2f kN", -f2(0)/1e3));
}

static void test_lumped_static()
{
    std::cout<<"lumped-mass line relaxed against the catenary"<<std::endl;
    model m;
    mooring_props p = chain();
    p.BA = 0.0;
    mooring_end an, fl;
    an.body = -1; an.p = Vec3(0.0, 0.0, -50.0);
    fl.body = -1; fl.p = Vec3(100.0, 0.0, -5.0);
    mooring_line line(m, an, fl, p, 40, nullptr, Vec3::Zero(), -50.0);
    line.initialize(m, 0.0);
    line.settle(m, 150.0);
    line.initialize(m, 0.0);
    Vec3 Fa, Fb;
    line.end_forces(Fa, Fb);
    catenary_line c;
    c.L = p.L; c.EA = p.EA; c.CB = 0.0;
    c.w = (p.mass - p.rho*0.25*M_PI*p.d*p.d)*9.81;
    double H = -1.0, V = 0.0;
    catenary_solve(c, 100.0, 45.0, H, V);
    // the seabed penalty lets the line sink by w / (kb d) = 1.4 mm; the end half-segment hangs on the fairlead
    check(std::fabs(-Fb(0) - H) < 0.01*H, fmt("fairlead horizontal force %.3f kN vs catenary %.3f kN", -Fb(0)/1e3, H/1e3));
    check(std::fabs(-Fb(2) - V) < 0.01*V, fmt("fairlead vertical force %.3f kN vs catenary %.3f kN", -Fb(2)/1e3, V/1e3));
    int on = 0;
    for(const Vec3 &r : line.nodes())
        if(r(2) < -50.0 + 0.01) ++on;
    const double LB = catenary_seabed_length(c, H, V);
    check(std::fabs(on*p.L/40.0 - LB) < 2.0*p.L/40.0, fmt("length on the seabed %.1f m vs catenary %.1f m", on*p.L/40.0, LB));
    // automatic time step: the smallest of the axial, seabed and damping limits of the explicit scheme
    const double mn = 50.0*3.0, dta = 0.5*std::sqrt(mn/(1.0e8/3.0)), dtb = 0.5*std::sqrt(mn/(3.0e6*0.1*3.0)), dtc = 0.5*mn/(3.0e5*0.1*3.0);
    check(std::fabs(line.dt_line() - std::min(dta, std::min(dtb, dtc))) < 1.0e-12, fmt("automatic line time step %.2e s", line.dt_line()));
}

static void test_surge_decay()
{
    std::cout<<"surge decay: dynamic line against quasi-static line"<<std::endl;
    double period[2], drift[2];
    for(int k=0; k<2; ++k)
    {
        model m;
        rigid_body b;
        b.mass = 2.0e5;
        b.inertia = Vec3(1e7, 1e7, 1e7).asDiagonal();
        m.add_body(b);
        m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 0, Vec3::Zero(), Vec3::UnitX()));
        mooring_props p = chain();
        p.BA = 5.0e4;
        p.cb = 3.0e4;
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(-100.0, 0.0, -50.0);
        fl.body = 0;  fl.p = Vec3(-5.0, 0.0, -5.0);
        mooring_end an2 = an, fl2 = fl;
        an2.p(0) = 100.0; fl2.p(0) = 5.0;
        if(k==0)
        {
            m.add_load(std::make_unique<mooring_qs>(m, an, fl, p));
            m.add_load(std::make_unique<mooring_qs>(m, an2, fl2, p));
            m.bodies[0].x(0) = 4.0;
        }
        else
        {
            // lines attached in the neutral position, then the body is offset and the lines relaxed there
            mooring_line &l1 = static_cast<mooring_line&>(m.add_load(std::make_unique<mooring_line>(m, an, fl, p, 12, nullptr, Vec3::Zero(), -50.0)));
            mooring_line &l2 = static_cast<mooring_line&>(m.add_load(std::make_unique<mooring_line>(m, an2, fl2, p, 12, nullptr, Vec3::Zero(), -50.0)));
            m.bodies[0].x(0) = 4.0;
            l1.settle(m, 100.0);
            l2.settle(m, 100.0);
        }
        timeint ti;
        ti.initialize(m, 0.0);
        double xprev = m.bodies[0].x(0);
        std::vector<double> zc;    // zero crossings (both directions)
        const double h = 0.2;
        for(int i=0; i<700; ++i)
        {
            ti.step(m, i*h, h);
            const double x = m.bodies[0].x(0);
            if(xprev*x < 0.0)
                zc.push_back(i*h + h*xprev/(xprev - x));
            xprev = x;
        }
        period[k] = zc.size() >= 3 ? 2.0*(zc.back() - zc.front())/(zc.size() - 1) : 0.0;
        drift[k] = ti.stats.failures;
    }
    check(period[0] > 0.0 && std::fabs(period[1] - period[0]) < 0.05*period[0],
          fmt("surge period: dynamic lines %.2f s, quasi-static %.2f s", period[1], period[0]));
    check(drift[0]==0 && drift[1]==0, "no failed steps");
}

static void test_drag()
{
    std::cout<<"line drag in still water"<<std::endl;
    // a taut horizontal line (no gravity) between two ground points, the inner nodes pushed sideways:
    // the drag on a node moving normal to the line is 1/2 rho Cdn d l |v| v
    model m;
    m.gravity = Vec3::Zero();
    mooring_props p = chain();
    p.L = 10.0;
    p.EA = 1.0e6;
    p.BA = 0.0;
    mooring_end a, b;
    a.body = -1; a.p = Vec3(0.0, 0.0, -10.0);
    b.body = -1; b.p = Vec3(10.0, 0.0, -10.0);
    mooring_line line(m, a, b, p, 10, nullptr, Vec3(0.0, 1.0, 0.0), -1.0e30, 1.0e-4);
    line.initialize(m, 0.0);
    // current of 1 m/s sideways: steady sideways force of the whole line = 1/2 rho Cdn d L U^2 once the
    // line has taken its deflected shape; check the total sideways end force after the transient
    line.settle(m, 30.0);
    line.initialize(m, 0.0);
    Vec3 Fa, Fb;
    line.end_forces(Fa, Fb);
    const double Fd = 0.5*p.rho*p.Cdn*p.d*p.L*1.0;
    // the line deflects, so the normal velocity is a bit less than U (tangential drag is small)
    check(std::fabs(Fa(1) + Fb(1) - Fd) < 0.05*Fd, fmt("sideways end forces %.1f N vs 1/2 rho Cd d L U^2 = %.1f N", Fa(1) + Fb(1), Fd));
}

static void test_shared_line()
{
    std::cout<<"line between two bodies (shared line): no seabed, hanging catenary"<<std::endl;
    // two bodies at the same height 80 m apart, 100 m of chain: no seabed by default, both lines hang
    // as the free catenary (before, the seabed at the lower end held the line up: no tension)
    model m;
    rigid_body b;
    b.mass = 1.0e5;
    b.inertia = Vec3(1e7, 1e7, 1e7).asDiagonal();
    b.x = Vec3(0.0, 0.0, -10.0);
    const int i0 = m.add_body(b);
    b.x = Vec3(80.0, 0.0, -10.0);
    const int i1 = m.add_body(b);
    mooring_props p = chain();
    p.L = 100.0;
    mooring_end a, e;
    a.body = i0; a.p = Vec3(0.0, 0.0, -10.0);
    e.body = i1; e.p = Vec3(80.0, 0.0, -10.0);
    check(default_seabed(a, e) < -1.0e29, "default seabed between two bodies: none");
    mooring_end g = a;
    g.body = -1;
    check(default_seabed(g, e)==-10.0, "default seabed of an anchored line: the anchor height");
    catenary_line c;
    c.L = p.L; c.EA = p.EA; c.CB = -1.0;
    c.w = (p.mass - p.rho*0.25*M_PI*p.d*p.d)*9.81;
    double H = -1.0, V = 0.0;
    catenary_solve(c, 80.0, 0.0, H, V);
    mooring_qs q(m, a, e, p);
    Vec3 Fa, Fb;
    q.forces(m, Fa, Fb);
    check(std::fabs(-Fb(0) - H) < 1.0e-6*H && std::fabs(Fa(0) - H) < 1.0e-6*H && std::fabs(-Fb(2) + Fa(2)) < 1.0e-6*H,
          fmt("quasi-static: H %.3f kN (catenary %.3f kN), symmetric", -Fb(0)/1e3, H/1e3));
    mooring_line l(m, a, e, p, 20, nullptr, Vec3::Zero(), default_seabed(a, e));
    l.initialize(m, 0.0);
    l.settle(m, 150.0);
    l.initialize(m, 0.0);
    l.end_forces(Fa, Fb);
    double zmin = 0.0;
    for(const Vec3 &r : l.nodes()) zmin = std::min(zmin, r(2));
    double xs, zs;
    catenary_profile(c, H, V, 0.5*p.L, xs, zs);
    check(std::fabs(-Fb(0) - H) < 0.01*H, fmt("dynamic: H %.3f kN vs catenary %.3f kN", -Fb(0)/1e3, H/1e3));
    check(std::fabs(zmin - (-10.0 + zs)) < 0.2, fmt("dynamic: sag %.2f m vs catenary %.2f m", -10.0 - zmin, -zs));
    // an anchor above the seabed (seabed given lower): the quasi-static line hangs free from the anchor
    mooring_end an, fl;
    an.body = -1; an.p = Vec3(0.0, 0.0, -40.0);
    fl.body = -1; fl.p = Vec3(60.0, 0.0, -5.0);
    mooring_props p2 = chain();
    p2.L = 75.0;
    mooring_qs q2(m, an, fl, p2, "raised", -50.0);
    q2.forces(m, Fa, Fb);
    catenary_line c2 = c;
    c2.L = p2.L;
    double H2 = -1.0, V2 = 0.0;
    catenary_solve(c2, 60.0, 35.0, H2, V2);
    check(std::fabs(-Fb(0) - H2) < 1.0e-6*H2 && std::fabs(-Fb(2) - V2) < 1.0e-6*V2,
          fmt("anchor 10 m above the seabed: hanging free, H %.3f kN (catenary %.3f kN)", -Fb(0)/1e3, H2/1e3));
}

static void test_wave_kinematics_per_step()
{
    std::cout<<"wave kinematics on the dynamic line: once per step against every Runge-Kutta stage"<<std::endl;
    // a chain between an anchor and a fixed fairlead 2 m below the surface in regular waves H 6 m, T 8 s:
    // the fairlead force with the kinematics at the converged nodes (start and end of the step, linear in
    // between) against the kinematics at every stage
    auto waves = std::make_shared<wave_field>(50.0);
    waves->add_regular(6.0, 8.0);
    waves->set_ramp(10.0);
    std::vector<double> F[2];
    double cpu[2];
    for(int k=0; k<2; ++k)
    {
        model m;
        mooring_props p = chain();
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(-100.0, 0.0, -50.0);
        fl.body = -1; fl.p = Vec3(0.0, 0.0, -2.0);
        mooring_line l(m, an, fl, p, 20, waves, Vec3::Zero(), -50.0);
        l.set_exact_kinematics(k==1);
        l.settle(m, 60.0);
        l.initialize(m, 0.0);
        const std::clock_t c0 = std::clock();
        const double h = 0.1;
        for(int i=1; i<=300; ++i)
        {
            l.commit(m, i*h, h);
            Vec3 Fa, Fb;
            l.end_forces(Fa, Fb);
            F[k].push_back(Fb.norm());
        }
        cpu[k] = double(std::clock() - c0)/CLOCKS_PER_SEC;
    }
    double dmax = 0.0, fmin = 1.0e300, fmax = -1.0e300;
    for(size_t i=0; i<F[0].size(); ++i)
    {
        dmax = std::max(dmax, std::fabs(F[0][i] - F[1][i]));
        fmin = std::min(fmin, F[1][i]);
        fmax = std::max(fmax, F[1][i]);
    }
    check(dmax < 0.01*(fmax - fmin), fmt("fairlead tension: largest difference %.3f kN of a range of %.2f kN", dmax/1e3, (fmax - fmin)/1e3));
    check(cpu[0] < 0.5*cpu[1], fmt("cost: %.3f s against %.3f s at every stage", cpu[0], cpu[1]));
}

static void test_end_motion_and_inertia()
{
    std::cout<<"dynamic line: Hermite end motion and the inertia of the end nodes"<<std::endl;
    // (a) a pretensioned line between two bodies accelerated along it (no gravity, no water): once the axial
    // ringing has decayed the line moves as a rigid rod and the bodies carry its whole mass,
    // F_a + F_b = -m L a (without the end nodes' inertia: -m (L - l) a)
    {
        model m;
        m.gravity = Vec3::Zero();
        rigid_body b;
        b.mass = 1.0;
        b.inertia = Mat3::Identity();
        const int i0 = m.add_body(b);
        b.x = Vec3(10.05, 0.0, 0.0);
        const int i1 = m.add_body(b);
        mooring_props p = chain();
        p.L = 10.0;
        p.EA = 1.0e6;
        p.BA = 2.0e4;
        p.rho = 0.0;
        mooring_end a, e;
        a.body = i0; a.p = Vec3(0.0, 0.0, 0.0);
        e.body = i1; e.p = Vec3(10.05, 0.0, 0.0);
        mooring_line l(m, a, e, p, 10, nullptr, Vec3::Zero(), -1.0e30);
        l.initialize(m, 0.0);
        const double acc = 0.2, h = 0.05;
        Vec3 Fa, Fb;
        for(int i=1; i<=400; ++i)
        {
            const double t = i*h;
            for(int k : {i0, i1})
            {
                m.bodies[k].x(0) = (k==i0 ? 0.0 : 10.05) + 0.5*acc*t*t;
                m.bodies[k].v(0) = acc*t;
                m.bodies[k].vd(0) = acc;
            }
            l.commit(m, t, h);
        }
        l.end_forces(Fa, Fb);
        const double F = Fa(0) + Fb(0), Fex = -p.mass*p.L*acc;
        check(std::fabs(F - Fex) < 0.005*std::fabs(Fex), fmt("rigid motion: F_a + F_b = %.3f N vs -m L a = %.3f N", F, Fex));
    }
    // (b) a chain on a fairlead heaving 2 m at T 8 s: the first harmonic of the fairlead tension does not depend
    // on the multibody step (the end velocity is continuous; linear end motion: 2.5 % at h 0.05 s)
    {
        double F1[2];
        const double hs[2] = {0.05, 0.0125};
        for(int k=0; k<2; ++k)
        {
            model m;
            rigid_body b;
            b.mass = 1.0;
            b.inertia = Mat3::Identity();
            const int ib = m.add_body(b);
            mooring_props p = chain();
            mooring_end an, fl;
            an.body = -1; an.p = Vec3(-100.0, 0.0, -50.0);
            fl.body = ib; fl.p = Vec3(0.0, 0.0, -5.0);
            mooring_line l(m, an, fl, p, 40, nullptr, Vec3::Zero(), -50.0);
            l.initialize(m, 0.0);
            l.settle(m, 150.0);
            l.initialize(m, 0.0);
            const double A = 2.0, w = 2.0*M_PI/8.0, h = hs[k];
            double sc = 0.0, ss = 0.0;
            int n = 0;
            const int N = (int)std::lround(64.0/h);
            for(int i=1; i<=N; ++i)
            {
                const double t = i*h;
                rigid_body &B = m.bodies[ib];
                B.x(2) = A*std::sin(w*t);
                B.v(2) = A*w*std::cos(w*t);
                B.vd(2) = -A*w*w*std::sin(w*t);
                l.commit(m, t, h);
                Vec3 Fa, Fb;
                l.end_forces(Fa, Fb);
                if(t > 32.0 + 1.0e-9)
                {
                    sc += Fb.norm()*std::cos(w*t);
                    ss += Fb.norm()*std::sin(w*t);
                    ++n;
                }
            }
            F1[k] = 2.0*std::hypot(sc, ss)/n;
        }
        check(std::fabs(F1[0] - F1[1]) < 0.002*F1[1], fmt("first harmonic of the fairlead tension %.3f kN (h 0.05 s), %.3f kN (h 0.0125 s)", F1[0]/1e3, F1[1]/1e3));
    }
}

static void test_surface_friction_guard()
{
    std::cout<<"dynamic line: free-surface crossing, seabed friction along and across, unstable line"<<std::endl;
    // (a) a vertical taut line from an anchor 20 m deep to a point 10 m above the water, 7 segments: the end
    // forces carry the weight in air less the buoyancy of exactly the submerged 20 m of line (counting
    // whole nodes in water or not: 19.3 m)
    {
        model m;
        mooring_props p = chain();
        p.L = 29.9;
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(0.0, 0.0, -20.0);
        fl.body = -1; fl.p = Vec3(0.0, 0.0, 10.0);
        mooring_line l(m, an, fl, p, 7, nullptr, Vec3::Zero(), -1.0e30);
        l.initialize(m, 0.0);
        l.settle(m, 30.0);
        l.initialize(m, 0.0);
        Vec3 Fa, Fb;
        l.end_forces(Fa, Fb);
        // the buoyancy of the unstretched length of 20 m of the stretched line (strain 30/29.9 - 1)
        const double A = 0.25*M_PI*p.d*p.d, Fex = -p.mass*p.L*9.81 + p.rho*A*20.0*(p.L/30.0)*9.81;
        check(std::fabs(Fa(2) + Fb(2) - Fex) < 0.002*p.rho*A*20.0*9.81,
              fmt("vertical end forces %.2f N vs -m L g + rho A 20 m g = %.2f N (buoyancy %.2f N)", Fa(2) + Fb(2), Fex, p.rho*A*20.0*9.81));
        check(std::fabs(l.tension(6) - std::fabs(Fb(2))) < 0.02*std::fabs(Fb(2)) && l.tension(6) > l.tension(0),
              fmt("end segment tensions Ta %.1f N < Tb %.1f N", l.tension(0), l.tension(6)));
    }
    // (b) a line on the seabed dragged by its ends at 0.5 m/s along and across itself: friction mu_axial
    // 0.5 along, mu 1.0 across, of the submerged weight (no drag); the ends at the depth the seabed penalty
    // lets the line sink (w / (kb d)). Stick-slip friction (default) and the regularized law (v_eps 0.01 m/s);
    // across, the pretensioned line rings sideways (+-25 %): the mean of the last 5 s
    double Fmean[2][2];
    for(int law=0; law<2; ++law)
    for(int dir=0; dir<2; ++dir)
    {
        mooring_props p = chain();
        p.L = 29.9;             // pretension 334 kN: the line bows little when dragged across
        p.Cdn = p.Cdt = 0.0;
        p.mu = 1.0;
        p.mu_axial = 0.5;
        p.v_eps = law==1 ? 0.01 : 0.0;
        const double w = (p.mass - p.rho*0.25*M_PI*p.d*p.d)*9.81, z = -50.0 - w/(p.kb*p.d);
        model m;
        rigid_body b;
        b.mass = 1.0;
        b.inertia = Mat3::Identity();
        b.x = Vec3(0.0, 0.0, z);
        const int i0 = m.add_body(b);
        b.x = Vec3(30.0, 0.0, z);
        const int i1 = m.add_body(b);
        mooring_end a, e;
        a.body = i0; a.p = Vec3(0.0, 0.0, z);
        e.body = i1; e.p = Vec3(30.0, 0.0, z);
        mooring_line l(m, a, e, p, 10, nullptr, Vec3::Zero(), -50.0);
        l.initialize(m, 0.0);
        l.settle(m, 20.0);
        l.initialize(m, 0.0);
        const double U = 0.5, h = 0.05;
        Vec3 Fa, Fb;
        double Fs = 0.0;
        for(int i=1; i<=400; ++i)
        {
            const double t = i*h;
            for(int k : {i0, i1})
            {
                m.bodies[k].x(dir) = (dir==0 && k==i1 ? 30.0 : 0.0) + U*t;
                m.bodies[k].v(dir) = U;
            }
            l.commit(m, t, h);
            l.end_forces(Fa, Fb);
            if(i > 300) Fs += (Fa(dir) + Fb(dir))/100.0;
        }
        Fmean[law][dir] = Fs;
        const double mu = dir==0 ? p.mu_axial : p.mu, Fex = -mu*w*p.L;
        const char *lw = law==0 ? "stick-slip" : "regularized";
        check(std::fabs(Fs - Fex) < (dir==0 ? 0.001 : 0.02)*std::fabs(Fex),
              std::string(lw) + fmt(dir==0 ? ": dragged along: friction %.1f N vs -mu_axial w L = %.1f N" : ": dragged across: friction %.1f N vs -mu w L = %.1f N",
                  Fs, Fex));
    }
    check(std::fabs(Fmean[0][1] - Fmean[1][1]) < 0.002*std::fabs(Fmean[1][1]),
          fmt("dragged across: stick-slip %.1f N, regularized %.1f N", Fmean[0][1], Fmean[1][1]));
    // (c) an explicit line step far too large: the run stops with the name of the line
    {
        model m;
        mooring_props p = chain();
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(-100.0, 0.0, -50.0);
        fl.body = -1; fl.p = Vec3(0.0, 0.0, -5.0);
        mooring_line l(m, an, fl, p, 20, nullptr, Vec3::Zero(), -50.0, 0.05, "unstable");
        l.initialize(m, 0.0);
        std::string msg;
        try
        {
            for(int i=1; i<=200; ++i)
                l.commit(m, 0.1*i, 0.1);
        }
        catch(const std::runtime_error &e) {msg = e.what();}
        check(msg.find("unstable") != std::string::npos && msg.find("dt_line") != std::string::npos, "unstable line: '" + msg + "'");
    }
}

static void test_static_start()
{
    std::cout<<"dynamic line: start in the static equilibrium of the segments"<<std::endl;
    // the chain of test_lumped_static without settling: from the static equilibrium of its 40 segments the
    // fairlead force is the catenary's and stays put with the ends held; from the continuous catenary
    // (no_static) the chords are shorter than the segments, the line drops and rings
    catenary_line c;
    mooring_props p = chain();
    c.L = p.L; c.EA = p.EA; c.CB = 0.0;
    c.w = (p.mass - p.rho*0.25*M_PI*p.d*p.d)*9.81;
    double H = -1.0, V = 0.0;
    catenary_solve(c, 100.0, 45.0, H, V);
    double dev[2];
    for(int k=0; k<2; ++k)
    {
        model m;
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(0.0, 0.0, -50.0);
        fl.body = -1; fl.p = Vec3(100.0, 0.0, -5.0);
        mooring_line l(m, an, fl, p, 40, nullptr, Vec3::Zero(), -50.0);
        l.set_static_start(k==0);
        l.initialize(m, 0.0);
        Vec3 Fa, Fb;
        l.end_forces(Fa, Fb);
        const double F0 = Fb.norm();
        if(k==0)
            check(std::fabs(-Fb(0) - H) < 0.005*H && std::fabs(-Fb(2) - V) < 0.005*V,
                  fmt("static start: fairlead force H %.3f kN, V %.3f kN (catenary %.3f, %.3f kN)", -Fb(0)/1e3, -Fb(2)/1e3, H/1e3, V/1e3));
        dev[k] = 0.0;
        for(int i=1; i<=200; ++i)
        {
            l.commit(m, 0.05*i, 0.05);
            l.end_forces(Fa, Fb);
            dev[k] = std::max(dev[k], std::fabs(Fb.norm() - F0)/F0);
        }
    }
    check(dev[0] < 0.002, fmt("ends held for 10 s: fairlead tension within %.2f %% of the start", 100.0*dev[0]));
    check(dev[1] > 10.0*dev[0], fmt("(from the continuous catenary: %.1f %%)", 100.0*dev[1]));
    // a stiff chain (EA 1.5e9 N, 140 kg/m): its catenary start has slack segments that the per-line Newton
    // cannot leave; the whole-system solve (compression stage, relaxation) finds the equilibrium
    {
        model m;
        mooring_props q = chain();
        q.d = 0.15; q.mass = 140.0; q.EA = 1.5e9; q.BA = 1.0e6;
        mooring_end an, fl;
        an.body = -1; an.p = Vec3(0.0, 0.0, -50.0);
        fl.body = -1; fl.p = Vec3(100.0, 0.0, -5.0);
        mooring_line l(m, an, fl, q, 40, nullptr, Vec3::Zero(), -50.0);
        const bool ok = l.equilibrium(m, 0.0);
        l.initialize(m, 0.0);
        Vec3 Fa, Fb;
        l.end_forces(Fa, Fb);
        catenary_line cs;
        cs.L = q.L; cs.EA = q.EA; cs.CB = 0.0;
        cs.w = (q.mass - q.rho*0.25*M_PI*q.d*q.d)*9.81;
        double Hs = -1.0, Vs = 0.0;
        catenary_solve(cs, 100.0, 45.0, Hs, Vs);
        check(ok && std::fabs(-Fb(0) - Hs) < 0.01*Hs, fmt("stiff chain: static start converged, H %.2f kN (catenary %.2f kN)", -Fb(0)/1e3, Hs/1e3));
    }
}

// sum of the converged forces of the lines on their body ends (ground included)
static Vec3 attachment_sum(const mooring_system &s)
{
    Vec3 F = Vec3::Zero();
    for(int k=0; k<s.nlines(); ++k)
    {
        Vec3 Fa, Fb;
        s.line_end_forces(k, Fa, Fb);
        // ends on points carry internal forces: only the body ends count (here all ends not on points)
        F += Fa + Fb;
    }
    return F;
}

static void test_points()
{
    std::cout<<"mooring systems: free points (clump weight, buoy, bridle)"<<std::endl;
    const double g = 9.81, rho = 1025.0;
    mooring_props p = chain();
    p.L = 60.0;
    const double w = (p.mass - rho*0.25*M_PI*p.d*p.d)*g;
    // (a) a clump weight hanging between two lines from two ground points 80 m apart (no seabed): in the static
    // equilibrium the ends carry the wet weight of the lines and the clump
    {
        model m;
        mooring_system s(m);
        mooring_point_props cp;
        cp.mass = 5000.0;
        cp.volume = 5000.0/7850.0;
        const int c = s.add_point("clump", Vec3(40.0, 0.0, -40.0), cp);
        mooring_end a, b, pc;
        a.p = Vec3(0.0, 0.0, -10.0);
        b.p = Vec3(80.0, 0.0, -10.0);
        pc.point = c;
        s.add_line(m, "left", a, pc, p, 20);
        s.add_line(m, "right", pc, b, p, 20);
        const bool ok = s.equilibrium(m, 0.0);
        s.initialize(m, 0.0);
        Vec3 Fa0, Fb0, Fa1, Fb1;
        s.line_end_forces(0, Fa0, Fb0);
        s.line_end_forces(1, Fa1, Fb1);
        const double Wc = (cp.mass - rho*cp.volume)*g, Wt = 2.0*w*p.L + Wc;
        check(ok && std::fabs(Fa0(2) + Fb1(2) + Wt) < 1.0e-4*Wt && std::fabs(Fa0(0) + Fb1(0)) < 1.0e-6*Wt,
              fmt("clump: ends carry %.2f kN of %.2f kN wet weight (lines %.2f, clump %.2f kN)", -(Fa0(2) + Fb1(2))/1e3, Wt/1e3, 2.0*w*p.L/1e3, Wc/1e3));
        check(std::fabs(s.point_position(c)(0) - 40.0) < 1.0e-6, fmt("clump in the middle: x %.6f m, z %.3f m", s.point_position(c)(0), s.point_position(c)(2)));
        // the point balances the two line ends and its own weight
        const Vec3 Fp = Fb0 + Fa1 + Vec3(0.0, 0.0, -Wc);
        check(Fp.norm() < 1.0e-6*Wt, fmt("forces on the clump balance to %.2e N", Fp.norm()));
    }
    // (b) a buoy (volume 4 m^3, 500 kg) on a riser from an anchor on the seabed, a second line on to a fixed
    // point: the ends carry the line weight less the buoy's net buoyancy, the riser lifts off the seabed
    {
        model m;
        mooring_system s(m);
        mooring_point_props bp;
        bp.mass = 500.0;
        bp.volume = 4.0;
        bp.Ca = 1.0;
        const int c = s.add_point("buoy", Vec3(10.0, 0.0, -25.0), bp);
        mooring_end a, b, pc;
        a.p = Vec3(0.0, 0.0, -50.0);
        b.p = Vec3(50.0, 0.0, -20.0);
        pc.point = c;
        mooring_props q = p;
        q.L = 30.0;
        s.add_line(m, "riser", a, pc, q, 15, -50.0);
        q.L = 42.0;
        s.add_line(m, "tail", pc, b, q, 15);
        const bool ok = s.equilibrium(m, 0.0);
        s.initialize(m, 0.0);
        Vec3 Fa0, Fb0, Fa1, Fb1;
        s.line_end_forces(0, Fa0, Fb0);
        s.line_end_forces(1, Fa1, Fb1);
        const double Wb = (bp.mass - rho*bp.volume)*g, Wt = w*72.0 + Wb;
        check(ok && std::fabs(Fa0(2) + Fb1(2) + Wt) < 1.0e-4*std::fabs(w*72.0) && Fa0(2) > 0.0,
              fmt("buoy: ends carry %.3f kN (lines %.3f kN, buoy %.3f kN), the riser lifts the anchor by %.3f kN", -(Fa0(2) + Fb1(2))/1e3, w*72.0/1e3, Wb/1e3, Fa0(2)/1e3));
    }
    // (c) bridle: a body held by two bridle lines to a point and a main line from the point to an anchor; the
    // forces on the body and the anchor balance the wet weight of the lines
    {
        model m;
        rigid_body b;
        b.mass = 1.0e5;
        b.inertia = Vec3(1e7, 1e7, 1e7).asDiagonal();
        const int ib = m.add_body(b);
        mooring_system s(m);
        const int c = s.add_point("bridle", Vec3(-15.0, 0.0, -10.0), mooring_point_props());
        mooring_end f1, f2, an, pc;
        f1.body = ib; f1.p = Vec3(-5.0, 4.0, -2.0);
        f2.body = ib; f2.p = Vec3(-5.0, -4.0, -2.0);
        an.p = Vec3(-120.0, 0.0, -50.0);
        pc.point = c;
        mooring_props q = p;
        q.L = 12.0;
        s.add_line(m, "bridle1", f1, pc, q, 6);
        s.add_line(m, "bridle2", f2, pc, q, 6);
        q.L = 115.0;
        s.add_line(m, "main", pc, an, q, 30, -50.0);
        const bool ok = s.equilibrium(m, 0.0);
        s.initialize(m, 0.0);
        Vec3 F = Vec3::Zero(), Fa, Fb;
        s.line_end_forces(0, Fa, Fb); F += Fa;
        s.line_end_forces(1, Fa, Fb); F += Fa;
        s.line_end_forces(2, Fa, Fb); F += Fb;
        // the seabed carries part of the main line: compare horizontally (no friction), and the bridle point
        check(ok && std::fabs(F(0)) < 1.0e-4*w*q.L && std::fabs(F(1)) < 1.0e-6*w*q.L,
              fmt("bridle: horizontal forces on body and anchor balance to %.2e N (body pulled by %.2f kN)", std::hypot(F(0), F(1)), 0.0));
        std::vector<int> bodies;
        s.coupled_bodies(bodies);
        check(bodies.size()==3 && s.npoints()==1 && std::fabs(s.point_position(c)(1)) < 1.0e-6,
              fmt("bridle point on the symmetry plane (y %.1e m), %.0f attachments", s.point_position(c)(1), (double)bodies.size()));
    }
    // (d) pendulum: a 1000 kg clump on a light 10 m line from a fixed point in air: period 2 pi sqrt(L/g)
    {
        model m;
        mooring_system s(m);
        mooring_point_props cp;
        cp.mass = 1000.0;
        mooring_props q = p;
        q.L = 10.0;
        q.mass = 1.0;
        q.d = 0.01;
        q.EA = 1.0e7;
        q.BA = 1.0e3;
        q.rho = 0.0;
        q.Cdn = q.Cdt = 0.0;
        const int c = s.add_point("bob", Vec3(1.0, 0.0, -std::sqrt(100.0 - 1.0)), cp);
        mooring_end top, pc;
        top.p = Vec3::Zero();
        pc.point = c;
        s.add_line(m, "string", top, pc, q, 10);
        s.set_static_start(false);
        s.initialize(m, 0.0);
        double xprev = s.point_position(c)(0);
        std::vector<double> zc;
        const double h = 0.01;
        for(int i=1; i<=1200; ++i)
        {
            s.commit(m, i*h, h);
            const double x = s.point_position(c)(0);
            if(xprev*x < 0.0)
                zc.push_back((i - 1)*h + h*xprev/(xprev - x));
            xprev = x;
        }
        const double T = zc.size() >= 3 ? 2.0*(zc.back() - zc.front())/(zc.size() - 1) : 0.0;
        // the pendulum length to the clump centre with the line's own mass: L_eff from the physical
        // pendulum of the clump and the uniform line (I/(m r)), small-amplitude correction (1 + th0^2/16)
        const double ml = q.mass*q.L, I = cp.mass*q.L*q.L + ml*q.L*q.L/3.0, mr = cp.mass*q.L + 0.5*ml*q.L;
        const double th0 = std::asin(0.1), Tex = 2.0*M_PI*std::sqrt(I/(mr*9.81))*(1.0 + th0*th0/16.0);
        check(std::fabs(T - Tex) < 0.003*Tex, fmt("pendulum period %.4f s vs %.4f s", T, Tex));
    }
    // (e) a clump weight resting on the seabed: contact force = wet weight, sunk by W / (kb pi r^2)
    {
        model m;
        mooring_system s(m);
        mooring_point_props cp;
        cp.mass = 8000.0;
        cp.volume = 1.0;
        const double rc = std::cbrt(0.75/M_PI);
        const int c = s.add_point("clump", Vec3(20.0, 0.0, -50.0 + rc), cp, -50.0);
        mooring_end a, pc;
        a.p = Vec3(0.0, 0.0, -50.0);
        pc.point = c;
        mooring_props q = p;
        q.L = 25.0;
        s.add_line(m, "slack", a, pc, q, 10, -50.0);
        const bool ok = s.equilibrium(m, 0.0);
        const double Wc = (cp.mass - rho*cp.volume)*g, dz = Wc/(cp.kb*M_PI*rc*rc);
        check(ok && std::fabs(s.point_position(c)(2) - (-50.0 + rc - dz)) < 0.05*dz,
              fmt("clump on the seabed: sunk %.4f m (W / (kb pi r^2) = %.4f m)", -50.0 + rc - s.point_position(c)(2), dz));
    }
    // (f) the junction of two lines (a point without volume) on the seabed: it has no contact of its own, the
    // line ends carry it like the line's nodes (sunk by w / (kb d)), and it does not shorten the sub-step
    {
        model m;
        mooring_system s(m), s1(m);
        const int c = s.add_point("junction", Vec3(20.0, 0.0, -50.0), mooring_point_props(), -50.0);
        mooring_end a, b, pc;
        a.p = Vec3(0.0, 0.0, -50.0);
        b.p = Vec3(40.0, 0.0, -50.0);
        pc.point = c;
        mooring_props q = p;
        q.L = 19.99;       // a little tension: the junction's place along the lines is determined
        s.add_line(m, "left", a, pc, q, 10, -50.0);
        s.add_line(m, "right", pc, b, q, 10, -50.0);
        mooring_end j;
        j.p = Vec3(20.0, 0.0, -50.0);
        s1.add_line(m, "left only", a, j, q, 10, -50.0);
        const bool ok = s.equilibrium(m, 0.0);
        const double dz = (q.mass - rho*0.25*M_PI*q.d*q.d)*g/(q.kb*q.d), zj = s.point_position(c)(2);
        check(ok, "junction on the seabed: static equilibrium converged");
        check(std::fabs(zj - (-50.0 - dz)) < 0.02*dz, fmt("junction on the seabed: sunk %.3e m (w / (kb d) = %.3e m)", -50.0 - zj, dz));
        check(s.dt_line() == s1.dt_line(), fmt("junction: sub-step %.3e s (of the lines alone %.3e s)", s.dt_line(), s1.dt_line()));
    }
}

static void test_break()
{
    std::cout<<"line failure"<<std::endl;
    // (a) a body held in surge by two opposite lines; the right line lets go of the body at t = 4 s: the body
    // drifts to the left anchor, the broken line no longer acts on it, its free end sinks
    auto build = [](model &m, std::vector<mooring_line*> &ls, const mooring_break *b)
    {
        rigid_body bd;
        bd.mass = 2.0e5;
        bd.inertia = Vec3(1e7, 1e7, 1e7).asDiagonal();
        const int ib = m.add_body(bd);
        m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, ib, Vec3::Zero(), Vec3::UnitX()));
        mooring_props p = chain();
        for(int k=0; k<2; ++k)
        {
            mooring_end an, fl;
            an.body = -1; an.p = Vec3(k==0 ? -100.0 : 100.0, 0.0, -50.0);
            fl.body = ib; fl.p = Vec3(k==0 ? -5.0 : 5.0, 0.0, -5.0);
            auto l = std::make_unique<mooring_line>(m, an, fl, p, 15, nullptr, Vec3::Zero(), -50.0, 0.0, k==0 ? "left" : "right");
            if(k==1 && b) l->set_break(0, *b);
            ls.push_back(&static_cast<mooring_line&>(m.add_load(std::move(l))));
        }
        return ib;
    };
    {
        model m;
        std::vector<mooring_line*> ls;
        mooring_break b;
        b.t = 4.0;
        const int ib = build(m, ls, &b);
        timeint ti;
        ti.initialize(m, 0.0);
        double x4 = 0.0, zfree4 = 0.0;
        bool steps = true;
        for(int i=0; i<200; ++i)
        {
            steps = steps && ti.step(m, i*0.05, 0.05);
            if(i==79) {x4 = m.bodies[ib].x(0); zfree4 = ls[1]->nodes().back()(2);}
        }
        std::vector<int> cb;
        ls[1]->coupled_bodies(cb);
        Vec3 Fa, Fb;
        ls[1]->end_forces(Fa, Fb);
        check(steps && ls[1]->line_broken(0) && std::fabs(ls[1]->line_break_time(0) - 4.0) < 1.0e-9 && cb.size()==1,
              fmt("broken at t = %.2f s, the line no longer on the body (%.0f attachment left)", ls[1]->line_break_time(0), (double)cb.size()));
        check(std::fabs(x4) < 1.0e-3 && m.bodies[ib].x(0) < -0.5 && ls[1]->nodes().back()(2) < zfree4 - 5.0,
              fmt("body at x %.3f m when it breaks, %.2f m 6 s later; the free end sank %.1f m", x4, m.bodies[ib].x(0), zfree4 - ls[1]->nodes().back()(2)));
        std::vector<std::string> ch;
        ls[1]->channel_names(ch);
        std::vector<double> cv;
        ls[1]->channel_values(m, 10.0, cv);
        check(ch.size()==10 && ch[9]=="right_broken" && cv.size()==10 && cv[9]==1.0, "channel right_broken = 1");
    }
    // (b) at a tension: the fairlead pulled away at 0.5 m/s, the line breaks at the first step at which the
    // tension of its end segment reaches 300 kN
    {
        model m;
        rigid_body bd;
        bd.mass = 1.0;
        bd.inertia = Mat3::Identity();
        const int ib = m.add_body(bd);
        mooring_props p = chain();
        mooring_end an, fl;
        an.p = Vec3(-100.0, 0.0, -50.0);
        fl.body = ib; fl.p = Vec3(0.0, 0.0, -5.0);
        mooring_line l(m, an, fl, p, 30, nullptr, Vec3::Zero(), -50.0);
        mooring_break b;
        b.T = 3.0e5;
        l.set_break(0, b);
        l.initialize(m, 0.0);
        double Tprev = 0.0, Tlast = 0.0, tb = -1.0;
        const double h = 0.05;
        for(int i=1; i<=1000 && tb < 0.0; ++i)
        {
            m.bodies[ib].x(0) = 0.5*i*h;
            m.bodies[ib].v(0) = 0.5;
            Tprev = Tlast;
            l.commit(m, i*h, h);
            Tlast = std::max(l.tension(0), l.tension(l.nseg() - 1));
            if(l.line_broken(0)) tb = i*h;
        }
        check(tb > 0.0 && Tprev < b.T && Tlast >= b.T && l.line_break_time(0)==tb,
              fmt("broke at t = %.2f s at %.1f kN (%.1f kN the step before)", tb, Tlast/1e3, Tprev/1e3));
    }
    // (c) restart after the failure: the checkpoint of the broken system continues bitwise as the run itself
    {
        model m1, m2;
        std::vector<mooring_line*> l1, l2;
        mooring_break b;
        b.t = 1.0;
        const int i1 = build(m1, l1, &b);
        build(m2, l2, &b);
        timeint t1;
        t1.initialize(m1, 0.0);
        for(int i=0; i<40; ++i) t1.step(m1, i*0.05, 0.05);
        std::vector<double> st;
        m1.save_state(st);
        std::vector<double> ts;
        t1.save_state(ts);
        timeint t2;
        t2.initialize(m2, 0.0);
        size_t k = 0, kt = 0;
        const bool ok = m2.load_state(st, k) && t2.load_state(m2, ts, kt);
        for(int i=40; i<80; ++i) {t1.step(m1, i*0.05, 0.05); t2.step(m2, i*0.05, 0.05);}
        check(ok && l2[1]->line_broken(0) && m2.bodies[i1].x(0)==m1.bodies[i1].x(0) && l2[1]->nodes().back()==l1[1]->nodes().back(),
              fmt("restart after the failure: body at %.6f m in both runs", m2.bodies[i1].x(0)));
    }
}

static void test_bathymetry()
{
    std::cout<<"seabed: plane, grid, the bed of a sampled fluid"<<std::endl;
    // (a) a clump weight on a frictionless 10 % slope, held up the slope by a light taut line parallel to it:
    // tension W sin(theta), contact force W cos(theta)
    {
        model m;
        mooring_system s(m);
        const double sl = 0.1, th = std::atan(sl);
        s.set_seabed(std::make_shared<seabed_plane>(-50.0, sl, 0.0));
        mooring_point_props cp;
        cp.mass = 2000.0;
        cp.radius = 0.3;
        const double x0 = 20.0;
        const Vec3 n(-std::sin(th), 0.0, std::cos(th));
        const Vec3 pc = Vec3(x0, 0.0, -50.0 + sl*x0) + 0.3*n;
        const int c = s.add_point("clump", pc, cp, -50.0);
        mooring_props q = chain();
        q.d = 0.02; q.mass = 0.5; q.EA = 1.0e7; q.BA = 1.0e4; q.L = 30.0;
        mooring_end top, e;
        top.p = pc + 30.2*Vec3(std::cos(th), 0.0, std::sin(th));
        e.point = c;
        s.add_line(m, "holder", e, top, q, 6, -1.0e30);
        const bool ok = s.equilibrium(m, 0.0);
        s.initialize(m, 0.0);
        Vec3 Fa, Fb;
        s.line_end_forces(0, Fa, Fb);
        const double W = cp.mass*9.81, T = -Fb.dot(Vec3(std::cos(th), 0.0, std::sin(th)));
        // the line's own weight (0.5 kg/m wet) adds a little to the tension at its top
        check(ok && std::fabs(T - W*std::sin(th)) < 0.01*W*std::sin(th),
              fmt("on a 10 %% slope: line tension %.1f N vs W sin(theta) = %.1f N", T, W*std::sin(th)));
    }
    // (b) a chain on a grid of the same plane gives the same line as on the plane itself
    {
        Vec3 F[2];
        for(int k=0; k<2; ++k)
        {
            model m;
            std::shared_ptr<seabed_surface> sb;
            if(k==0) sb = std::make_shared<seabed_plane>(-50.0, 0.05, -0.02);
            else
            {
                std::vector<double> z;
                for(int j=0; j<11; ++j) for(int i=0; i<21; ++i) z.push_back(-50.0 + 0.05*(-150.0 + 10.0*i) - 0.02*(-50.0 + 10.0*j));
                sb = std::make_shared<seabed_grid>(21, 11, -150.0, -50.0, 10.0, 10.0, z);
            }
            mooring_props p = chain();
            mooring_end an, fl;
            an.p = Vec3(-100.0, 0.0, -50.0 - 5.0);
            fl.p = Vec3(0.0, 0.0, -5.0);
            mooring_line l(m, an, fl, p, 30, nullptr, Vec3::Zero(), -55.0);
            l.set_seabed(sb);
            l.initialize(m, 0.0);
            Vec3 Fa;
            l.end_forces(Fa, F[k]);
        }
        check((F[0] - F[1]).norm() < 1.0e-6*F[0].norm(), fmt("grid of a plane = plane: fairlead force %.4f kN, difference %.1e N", F[0].norm()/1e3, (F[0] - F[1]).norm()));
    }
    // (c) the bed of a resolved fluid at the nodes (sampled kinematics): 10 m above the line's flat seabed
    {
        model m;
        auto sk = std::make_shared<sampled_kinematics>();
        mooring_props p = chain();
        mooring_end an, fl;
        an.p = Vec3(-100.0, 0.0, -40.0);
        fl.p = Vec3(0.0, 0.0, -5.0);
        auto &l = static_cast<mooring_line&>(m.add_load(std::make_unique<mooring_line>(m, an, fl, p, 30, sk, Vec3::Zero(), -50.0)));
        l.set_seabed(nullptr, true);
        std::vector<Vec3> pts;
        sk->collect(m, pts);
        const int n = (int)pts.size();
        sk->add_sample(0.0, std::vector<Vec3>(n, Vec3::Zero()), std::vector<double>(n, 0.0), std::vector<int>(n, 1), std::vector<double>(n, -40.0));
        l.initialize(m, 0.0);
        double zmin = 0.0;
        for(const Vec3 &r : l.nodes()) zmin = std::min(zmin, r(2));
        check(zmin > -40.01 && zmin < -39.99, fmt("on the fluid's bed at -40 m (flat seabed -50 m): lowest node %.4f m", zmin));
    }
}

static void test_bending()
{
    std::cout<<"bending stiffness: cantilever, taut line"<<std::endl;
    // a 10 m cantilever in air (10 kg/m, EI 1e6 N m^2) clamped horizontally at x = 0, its free end on a point
    // without mass: tip deflection w L^4 / (8 EI), moment on the clamping body w L^2 / 2
    const double Lb = 10.0, mb = 10.0, EI = 1.0e6, w = mb*9.81;
    auto build = [&](model &m, mooring_system &s, int body, bool stat)
    {
        mooring_props q = chain();
        q.L = Lb; q.d = 0.05; q.mass = mb; q.EA = 1.0e9; q.BA = 0.0; q.rho = 0.0; q.EI = EI;
        q.Cdn = q.Cdt = 0.0;
        const int tip = s.add_point("tip", Vec3(Lb, 0.0, 0.0), mooring_point_props());
        mooring_end a, b;
        a.body = body; a.p = Vec3::Zero();
        b.point = tip;
        s.add_line(m, "beam", a, b, q, 20);
        s.set_clamp(0, 0, m, Vec3::UnitX());
        s.set_static_start(stat);
        return tip;
    };
    {
        model m;
        mooring_system s(m);
        const int tip = build(m, s, -1, true);
        const bool ok = s.equilibrium(m, 0.0);
        const double d = -s.point_position(tip)(2), dex = w*std::pow(Lb, 4)/(8.0*EI);
        check(ok && std::fabs(d - dex) < 0.01*dex, fmt("tip deflection %.5f m vs w L^4/(8 EI) = %.5f m", d, dex));
    }
    // a taut line with bending: first natural frequency f_string sqrt(1 + pi^2 EI/(T L^2)) (pinned ends;
    // a tension-only line cannot carry a free cantilever's compression in motion); excited by a step of one
    // end, the frequency of the largest response of the middle node
    for(int k=0; k<2; ++k)
    {
        model m;
        m.gravity = Vec3::Zero();
        rigid_body bd;
        bd.mass = 1.0;
        bd.inertia = Mat3::Identity();
        const int ib = m.add_body(bd);
        mooring_system s(m);
        mooring_props q = chain();
        q.d = 0.05; q.mass = 10.0; q.EA = 1.0e7; q.BA = 0.0; q.rho = 0.0; q.EI = k==0 ? 0.0 : 1.0e4;
        q.Cdn = q.Cdt = 0.0;
        q.L = 10.0/(1.0 + 1.0e4/q.EA);            // pretension 1e4 N over the span of 10 m
        mooring_end a, b;
        a.p = Vec3::Zero();
        b.body = ib; b.p = Vec3(10.0, 0.0, 0.0);
        s.add_line(m, "string", a, b, q, 20);
        s.initialize(m, 0.0);
        const double h = 0.005;
        std::vector<double> z;
        for(int i=1; i<=4000; ++i)
        {
            m.bodies[ib].x(2) = 0.002;
            s.commit(m, i*h, h);
            z.push_back(s.line_nodes(0)[10](2));
        }
        double zm = 0.0;
        for(double x : z) zm += x/z.size();
        for(double &x : z) x -= zm;
        double fbest = 0.0, pbest = -1.0;
        for(double f=1.2; f<=2.2; f+=0.0005)
        {
            double c = 0.0, sn = 0.0;
            for(size_t i=0; i<z.size(); ++i) {c += z[i]*std::cos(2.0*M_PI*f*i*h); sn += z[i]*std::sin(2.0*M_PI*f*i*h);}
            const double pw = c*c + sn*sn;
            if(pw > pbest) {pbest = pw; fbest = f;}
        }
        // the 20 lumped masses: (1/pi) sqrt(T/(m l^2)) sin(pi/(2 N)) for the string (0.09 % below 1/(2L) sqrt(T/m))
        const double T = 1.0e4, fs = std::sqrt(T/(q.mass*0.25))*std::sin(M_PI/40.0)/M_PI;
        const double fex = fs*std::sqrt(1.0 + M_PI*M_PI*q.EI/(T*100.0));
        check(std::fabs(fbest - fex) < 0.003*fex, fmt(k==0 ? "taut string: first frequency %.4f Hz vs %.4f Hz" : "taut line with EI 1e4 N m^2: first frequency %.4f Hz vs %.4f Hz", fbest, fex));
    }
    {
        model m;
        rigid_body bd;
        bd.mass = 1.0e3;
        bd.inertia = Mat3::Identity()*1.0e3;
        const int ib = m.add_body(bd);
        mooring_system s(m);
        build(m, s, ib, true);
        s.initialize(m, 0.0);
        VecX f = VecX::Zero(m.ndof());
        s.eval(m, 0.0, f);
        const double My = f(4), Mex = w*Lb*Lb/2.0;
        check(std::fabs(My - Mex) < 0.01*Mex && std::fabs(f(2) + w*Lb) < 0.01*w*Lb,
              fmt("clamping body: moment %.1f N m vs w L^2/2 = %.1f N m, force %.1f N", My, Mex, f(2)));
    }
}

static void test_damping_law()
{
    std::cout<<"axial damping law: taut segments only, or in every segment (BA_slack)"<<std::endl;
    // (a) a slack straight line of two segments, no gravity, no water, its body end moving at 1 m/s towards the
    // ground end: with BA_slack the end segment pushes the body away by BA v / l, else nothing
    for(int law=0; law<2; ++law)
    {
        model m;
        m.gravity = Vec3::Zero();
        rigid_body b;
        b.mass = 1.0;
        b.inertia = Mat3::Identity();
        b.x = Vec3(9.0, 0.0, 0.0);
        b.v = Vec3(-1.0, 0.0, 0.0);
        const int ib = m.add_body(b);
        mooring_props p = chain();
        p.L = 10.0;
        p.EA = 1.0e6;
        p.BA = 2.0e4;
        p.rho = 0.0;
        p.BA_slack = law==1;
        mooring_end a, e;
        a.p = Vec3::Zero();
        e.body = ib; e.p = b.x;
        mooring_line l(m, a, e, p, 2, nullptr, Vec3::Zero(), -1.0e30);
        l.set_static_start(false);
        l.initialize(m, 0.0);
        Vec3 Fa, Fb;
        l.end_forces(Fa, Fb);
        const double Fex = law==1 ? p.BA*1.0/(p.L/2) : 0.0;
        check(std::fabs(Fb(0) - Fex) < 1.0e-9*p.BA && std::fabs(l.tension(1) + Fex) < 1.0e-9*p.BA,
              fmt(law==1 ? "BA_slack: a compressed slack segment pushes the body by %.1f N, tension %.1f N (BA v / l = %.1f N)"
                         : "default: a compressed slack segment, force on the body %.1f N, tension %.1f N (%.1f N)",
                  Fb(0), l.tension(1), Fex));
    }
    // (b) a taut chain without a seabed on a fairlead heaving 0.5 m at T 8 s, all segments taut all the time:
    // both laws give the same line bit for bit
    std::vector<double> Tb[2];
    double Tmin = 1.0e300;
    for(int law=0; law<2; ++law)
    {
        model m;
        rigid_body b;
        b.mass = 1.0;
        b.inertia = Mat3::Identity();
        b.x = Vec3(0.0, 0.0, -5.0);
        const int ib = m.add_body(b);
        mooring_props p = chain();
        p.L = 109.0;
        p.BA = 2.0e5;
        p.BA_slack = law==1;
        mooring_end a, e;
        a.p = Vec3(-100.0, 0.0, -50.0);
        e.body = ib; e.p = b.x;
        mooring_line l(m, a, e, p, 20, nullptr, Vec3::Zero(), -1.0e30);
        l.initialize(m, 0.0);
        const double h = 0.05, w = 2.0*M_PI/8.0;
        for(int i=1; i<=320; ++i)
        {
            const double t = i*h;
            m.bodies[ib].x(2) = -5.0 + 0.5*std::sin(w*t);
            m.bodies[ib].v(2) = 0.5*w*std::cos(w*t);
            m.bodies[ib].vd(2) = -0.5*w*w*std::sin(w*t);
            l.commit(m, t, h);
            Tb[law].push_back(l.tension(l.nseg() - 1));
            if(law==0)
                for(int j=0; j<l.nseg(); ++j) Tmin = std::min(Tmin, l.tension(j));
        }
    }
    double dmax = 0.0;
    for(size_t i=0; i<Tb[0].size(); ++i) dmax = std::max(dmax, std::fabs(Tb[0][i] - Tb[1][i]));
    check(Tmin > 0.0 && dmax == 0.0, fmt("taut line: the laws differ by %.1e N (smallest segment tension %.1f N)", dmax, Tmin));
}

static void test_stick_slip()
{
    std::cout<<"seabed friction: stick-slip"<<std::endl;
    // (a) a chain lying straight along x on a seabed sloping 1:5 across it (z = -50 + 0.2 y), its ends held by two
    // bodies: with mu 0.5 > tan(theta) it stays where it lies (stick-slip); the regularized law (v_eps 0.01 m/s)
    // creeps downhill under the slope component of its weight; without friction it slides down to its sag
    auto slope = [](double vep, double mu, double T, mooring_system **keep, model *mk) -> double
    {
        model local;
        model &m = mk ? *mk : local;
        rigid_body b;
        b.mass = 1.0;
        b.inertia = Mat3::Identity();
        b.x = Vec3(0.0, 0.0, -50.0);
        const int i0 = m.add_body(b);
        b.x = Vec3(40.0, 0.0, -50.0);
        const int i1 = m.add_body(b);
        mooring_props p = chain();
        p.L = 40.0;
        p.mu = mu;
        p.v_eps = vep;
        auto *s = new mooring_system(m);
        s->set_seabed(std::make_shared<seabed_plane>(-50.0, 0.0, 0.2));
        mooring_end a, e;
        a.body = i0; a.p = Vec3(0.0, 0.0, -50.0);
        e.body = i1; e.p = Vec3(40.0, 0.0, -50.0);
        s->add_line(m, "on the slope", a, e, p, 20, -50.0);
        s->set_static_start(false);      // straight along x where it was laid (the frictionless equilibrium sags)
        s->initialize(m, 0.0);
        const double h = 0.05;
        for(int i=1; i<=(int)std::lround(T/h); ++i)
            s->commit(m, i*h, h);
        const double y = s->line_nodes(0)[10](1);
        if(keep) *keep = s; else delete s;
        return y;
    };
    const double ys = slope(0.0, 0.5, 30.0, nullptr, nullptr), yr = slope(0.01, 0.5, 30.0, nullptr, nullptr),
                 y0 = slope(0.0, 0.0, 30.0, nullptr, nullptr);
    check(std::fabs(ys) < 0.005 && yr < -0.05 && y0 < yr,
          fmt("side slope 1:5, mu 0.5: the middle moves %.4f m in 30 s (regularized %.3f m, frictionless %.3f m)", ys, yr, y0));
    // (b) a clump weight (2 t, 0.25 m^3, mu 0.6) on a flat seabed pulled by a neutrally buoyant rope from a body
    // moving at 0.05 m/s: it breaks away at mu_static mu W (1.5) and slides at mu W
    for(int k=0; k<2; ++k)
    {
        const double mus = k==0 ? 1.5 : 1.0;
        model m;
        rigid_body b;
        b.mass = 1.0;
        b.inertia = Mat3::Identity();
        mooring_point_props cp;
        cp.mass = 2000.0;
        cp.volume = 0.25;
        cp.mu = 0.6;
        cp.mu_static = mus;
        const double rc = std::cbrt(0.75*cp.volume/M_PI), zc = -50.0 + rc;
        b.x = Vec3(20.0, 0.0, zc);
        const int ib = m.add_body(b);
        mooring_system s(m);
        const int c = s.add_point("clump", Vec3(0.0, 0.0, zc), cp, -50.0);
        mooring_props q;
        q.d = 0.05;
        q.mass = 1025.0*0.25*M_PI*q.d*q.d;      // neutral: straight
        q.EA = 1.0e6;
        q.BA = 1.0e3;
        q.Cdn = q.Cdt = q.Can = 0.0;
        q.L = 20.0;
        mooring_end pc, e;
        pc.point = c;
        e.body = ib; e.p = b.x;
        s.add_line(m, "rope", pc, e, q, 10);
        s.initialize(m, 0.0);
        const double h = 0.05, U = 0.05, W = (cp.mass - 1025.0*cp.volume)*9.81;
        double Tmax = 0.0, Tsl = 0.0;
        int nsl = 0;
        for(int i=1; i<=1200; ++i)
        {
            const double t = i*h;
            m.bodies[ib].x(0) = 20.0 + U*t;
            m.bodies[ib].v(0) = U;
            s.commit(m, t, h);
            const double T = s.line_tension(0, 0);
            Tmax = std::max(Tmax, T);
            if(t > 40.0) {Tsl += T; ++nsl;}
        }
        Tsl /= nsl;
        if(k==0)
            check(std::fabs(Tmax - mus*cp.mu*W) < 0.03*mus*cp.mu*W && s.point_position(c)(0) > 0.5,
                  fmt("clump breaks away at %.0f N (mu_static mu W = %.0f N), slides %.2f m", Tmax, mus*cp.mu*W, s.point_position(c)(0)));
        else
            check(std::fabs(Tsl - cp.mu*W) < 0.01*cp.mu*W,
                  fmt("clump sliding: rope tension %.0f N (mu W = %.0f N)", Tsl, cp.mu*W));
    }
    // (c) restart: the side-slope chain with one end swept 0.5 m across at T 6 s (sticking and sliding nodes),
    // saved at 10 s: the restarted run continues bit for bit
    {
        auto run = [&](bool restart, std::vector<double> *mid) -> Vec3
        {
            model m;
            mooring_system *s = nullptr;
            slope(0.0, 0.5, 0.0, &s, &m);
            const double h = 0.05;
            std::vector<double> st;
            for(int i=1; i<=400; ++i)
            {
                const double t = i*h;
                m.bodies[1].x(1) = 0.5*std::sin(2.0*M_PI*t/6.0);
                m.bodies[1].v(1) = 0.5*2.0*M_PI/6.0*std::cos(2.0*M_PI*t/6.0);
                s->commit(m, t, h);
                if(i==200)
                {
                    if(!restart) {s->save_state(st); *mid = st;}
                    else
                    {
                        // a fresh system from the checkpoint
                        model m2;
                        mooring_system *s2 = nullptr;
                        slope(0.0, 0.5, 0.0, &s2, &m2);
                        size_t k = 0;
                        s2->load_state(*mid, k);
                        delete s;
                        s = s2;
                        m.bodies[1].x(1) = 0.5*std::sin(2.0*M_PI*t/6.0);
                    }
                }
            }
            const Vec3 r = s->line_nodes(0)[10];
            delete s;
            return r;
        };
        std::vector<double> mid;
        const Vec3 r1 = run(false, &mid), r2 = run(true, &mid);
        check(r1==r2 && mid.size() > 0, fmt("restart with stick-slip anchors: the middle node at y %.6f m in both runs (%.0f values saved)",
                                            r1(1), (double)mid.size()));
    }
}

int main()
{
    test_catenary();
    test_quasi_static();
    test_lumped_static();
    test_surge_decay();
    test_drag();
    test_shared_line();
    test_wave_kinematics_per_step();
    test_end_motion_and_inertia();
    test_surface_friction_guard();
    test_static_start();
    test_points();
    test_break();
    test_bathymetry();
    test_bending();
    test_damping_law();
    test_stick_slip();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
