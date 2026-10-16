// Standalone tests of the REEF3D::MAST mooring lines (src/mast/mooring): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         mast_mooring_test.cpp ../../src/mast/core/*.cpp ../../src/mast/structure/*.cpp
//         ../../src/mast/hydro/*.cpp ../../src/mast/mooring/*.cpp -o mast_mooring_test
// Run:    ./mast_mooring_test
//
// Tests: elastic catenary (fairlead forces from the end position and back, hanging free and with a
// seabed portion and friction, slack line with length to spare, profile ending at the fairlead,
// anchor forces), quasi-static line on a
// body (force direction, symmetric spread), lumped-mass line relaxed between fixed ends against the
// catenary (fairlead tension, seabed length), surge decay of a moored body with the dynamic line
// against the quasi-static line (period), drag of a taut line in a current.
// Architect: Hans Bihs
#include"mast/core/mast_model.h"
#include"mast/core/mast_timeint.h"
#include"mast/core/mast_gravity.h"
#include"mast/structure/mast_linear6.h"
#include"mast/mooring/mast_catenary.h"
#include"mast/mooring/mast_mooring.h"
#include<iostream>
#include<string>
#include<cstdio>
#include<cmath>
#include<memory>

using namespace mast;

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

int main()
{
    test_catenary();
    test_quasi_static();
    test_lumped_static();
    test_surge_decay();
    test_drag();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
