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
// the stop of an unstable line, the start from the static equilibrium of the segments.
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_gravity.h"
#include"coral/structure/coral_linear6.h"
#include"coral/mooring/coral_catenary.h"
#include"coral/mooring/coral_mooring.h"
#include"coral/hydro/coral_waves.h"
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
    // lets the line sink (w / (kb d))
    for(int dir=0; dir<2; ++dir)
    {
        mooring_props p = chain();
        p.L = 29.9;             // pretension 334 kN: the line bows little when dragged across
        p.Cdn = p.Cdt = 0.0;
        p.mu = 1.0;
        p.mu_axial = 0.5;
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
        for(int i=1; i<=400; ++i)
        {
            const double t = i*h;
            for(int k : {i0, i1})
            {
                m.bodies[k].x(dir) = (dir==0 && k==i1 ? 30.0 : 0.0) + U*t;
                m.bodies[k].v(dir) = U;
            }
            l.commit(m, t, h);
        }
        l.end_forces(Fa, Fb);
        const double mu = dir==0 ? p.mu_axial : p.mu, Fex = -mu*w*p.L;
        check(std::fabs(Fa(dir) + Fb(dir) - Fex) < 0.01*std::fabs(Fex),
              fmt(dir==0 ? "dragged along: friction %.1f N vs -mu_axial w L = %.1f N" : "dragged across: friction %.1f N vs -mu w L = %.1f N",
                  Fa(dir) + Fb(dir), Fex));
    }
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
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
