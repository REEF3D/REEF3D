// Standalone tests of the REEF3D::CORAL co-rotational beams (src/coral/structure, core static solver):
// no MPI, no REEF3D.
// Build:  see coral_tests.sh (links all CORAL objects)
// Run:    ./coral_beam_test
//
// Tests: cantilever under a tip load (Timoshenko and Euler-Bernoulli), axial and torsion; roll-up of a
// cantilever into a circle under the tip moment 2 pi EI/L (large rotations, load steps); Bathe's 45°
// bend under an out-of-plane tip load against the published tip position; first bending frequency of
// a cantilever with lumped mass; free beam in rigid motion (no strain, energy kept); a beam end on a
// rigid body at an offset (moment from the offset).
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_static.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_constraint.h"
#include"coral/structure/coral_beam.h"
#include"coral/structure/coral_frame.h"
#include<cmath>
#include<cstdio>
#include<functional>
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

// force F (global) and moment M (global) on body b, set from outside
class point_load : public load
{
public:
    point_load(int b, const Vec3 &F, const Vec3 &M) : b_(b), F(F), M(M) {}
    std::string name() const override {return "point_load";}
    coupling_kind coupling() const override {return coupling_kind::local;}
    void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
    void eval(const model &m, double, VecX &f) const override
    {
        m.add_force(f, b_, m.body(b_).x, F);
        m.add_moment(f, b_, M);
    }
    int b_;
    Vec3 F, M;
};

static const double E = 2.1e11, G = 8.08e10, rho = 7850.0;

// cantilever along x from the ground at the origin, n elements; returns the tip node
static int cantilever(model &m, double L, int n, const beam_section &s, std::vector<int> *st=nullptr)
{
    frame_builder fb(m);
    member_spec ms;
    ms.name = "c";
    ms.p1 = Vec3::Zero();
    ms.p2 = Vec3(L, 0.0, 0.0);
    ms.n = n;
    ms.s1 = s;
    ms.body1 = -1;
    const std::vector<int> v = fb.add_member(ms);
    fb.finish();
    if(st) *st = v;
    return v.back();
}

static void test_cantilever()
{
    std::cout<<"cantilever: tip load, axial, torsion"<<std::endl;
    const double L = 10.0;
    const beam_section s = beam_section::tube(0.5, 0.02, E, G, rho);
    // transverse tip load (small): PL^3/(3EI) + PL/(GA)
    for(int timo=1; timo>=0; --timo)
    {
        beam_section ss = s;
        if(!timo) ss.GAy = ss.GAz = 0.0;
        model m;
        m.gravity = Vec3::Zero();
        const int tip = cantilever(m, L, 10, ss);
        const double P = 100.0;
        m.add_load(std::make_unique<point_load>(tip, Vec3(0, 0, -P), Vec3::Zero()));
        const static_result r = static_solve(m, 0.0);
        const double w = -m.bodies[tip].x(2);
        const double wex = P*L*L*L/(3.0*ss.EIy) + (timo ? P*L/ss.GAz : 0.0);
        check(r.converged && std::fabs(w - wex) < 1.0e-5*wex, fmt(timo ? "Timoshenko tip deflection %.6e m vs %.6e m" : "Euler-Bernoulli tip deflection %.6e m vs %.6e m", w, wex));
        // tip rotation PL^2/(2EI)
        const double th = log_so3(m.bodies[tip].R)(1);
        check(std::fabs(th - P*L*L/(2.0*ss.EIy)) < 1.0e-5*P*L*L/(2.0*ss.EIy), fmt("tip rotation %.6e rad vs PL^2/(2EI) %.6e rad", th, P*L*L/(2.0*ss.EIy)));
    }
    {
        model m;
        m.gravity = Vec3::Zero();
        const int tip = cantilever(m, L, 4, s);
        const double P = 1.0e5, T = 2.0e4;
        m.add_load(std::make_unique<point_load>(tip, Vec3(P, 0, 0), Vec3(T, 0, 0)));
        static_solve(m, 0.0);
        const double u = m.bodies[tip].x(0) - L;
        const double phi = log_so3(m.bodies[tip].R)(0);
        check(std::fabs(u - P*L/s.EA) < 1.0e-6*P*L/s.EA, fmt("axial: %.6e m vs PL/EA %.6e m", u, P*L/s.EA));
        check(std::fabs(phi - T*L/s.GJ) < 1.0e-6*T*L/s.GJ, fmt("torsion: %.6e rad vs TL/GJ %.6e rad", phi, T*L/s.GJ));
    }
}

static void test_rollup()
{
    std::cout<<"roll-up under a tip moment (large rotations)"<<std::endl;
    const double L = 10.0;
    beam_section s;
    s.EA = 1.0e8; s.EIy = s.EIz = 1.0e4; s.GJ = 1.0e4; s.m = 1.0; s.rIx = s.rIy = s.rIz = 1.0e-3;
    model m;
    m.gravity = Vec3::Zero();
    const int tip = cantilever(m, L, 20, s);
    point_load &pl = static_cast<point_load&>(m.add_load(std::make_unique<point_load>(tip, Vec3::Zero(), Vec3::Zero())));
    const double Mmax = 2.0*M_PI*s.EIy/L;
    bool ok = true;
    for(int k=1; k<=20; ++k)
    {
        pl.M = Vec3(0.0, -Mmax*k/20.0, 0.0);
        ok = ok && static_solve(m, 0.0).converged;
        if(k==10)
        {
            // half the moment: a half circle, tip at (0, 0, 2 L/pi) ... radius EI/M = L/pi
            const Vec3 tip_ex(0.0, 0.0, 2.0*L/M_PI);
            check((m.bodies[tip].x - tip_ex).norm() < 2.0e-3*L, fmt("half moment: tip at (%.4f, %.4f) vs (0, %.4f) (half circle)", m.bodies[tip].x(0), m.bodies[tip].x(2), 2.0*L/M_PI));
        }
    }
    check(ok, "all load steps converged");
    check(m.bodies[tip].x.norm() < 2.0e-3*L, fmt("full moment 2 pi EI/L: tip back at the root (|x| = %.2e m)", m.bodies[tip].x.norm()));
    check((m.bodies[tip].R - Mat3::Identity()).norm() < 1.0e-2, fmt("tip turned by 2 pi (|R - I| = %.1e)", (m.bodies[tip].R - Mat3::Identity()).norm()));
}

static void test_bend45()
{
    std::cout<<"45 degree bend under an out-of-plane tip load (Bathe & Bolourchi)"<<std::endl;
    // radius 100 in the x-y plane, centre (0, 100, 0), from the origin along +x; square section 1 x 1,
    // E = 1e7, G = E/2; 8 straight elements; tip load 600 in +z
    beam_section s;
    s.EA = 1.0e7; s.EIy = s.EIz = 1.0e7/12.0; s.GJ = 0.5e7*0.1406; s.m = 1.0; s.rIx = s.rIy = s.rIz = 1.0e-2;
    model m;
    m.gravity = Vec3::Zero();
    frame_builder fb(m);
    const int n = 8;
    std::vector<int> nodes;
    int prev = -1;
    Vec3 pprev = Vec3::Zero();
    for(int i=1; i<=n; ++i)
    {
        const double th = M_PI/4.0*i/n;
        const Vec3 p(100.0*std::sin(th), 100.0 - 100.0*std::cos(th), 0.0);
        member_spec ms;
        ms.name = "b" + std::to_string(i);
        ms.p1 = pprev; ms.p2 = p; ms.n = 1; ms.s1 = s;
        ms.orient = Vec3::UnitZ();
        if(i==1) ms.body1 = -1;
        const std::vector<int> st = fb.add_member(ms);
        prev = st.back();
        pprev = p;
    }
    fb.finish();
    point_load &pl = static_cast<point_load&>(m.add_load(std::make_unique<point_load>(prev, Vec3::Zero(), Vec3::Zero())));
    bool ok = true;
    for(int k=1; k<=20; ++k)
    {
        pl.F = Vec3(0.0, 0.0, 600.0*k/20.0);
        ok = ok && static_solve(m, 0.0).converged;
    }
    const Vec3 x = m.bodies[prev].x, ref(47.23, 15.79, 53.37);
    check(ok, "load steps converged");
    check((x - ref).norm() < 0.01*ref.norm(), fmt("tip (%.2f, %.2f, %.2f) vs Simo & Vu-Quoc (47.23, 15.79, 53.37)", x(0), x(1), x(2)));
}

static void test_frequency()
{
    std::cout<<"cantilever bending frequency (lumped mass)"<<std::endl;
    const double L = 10.0;
    beam_section s = beam_section::tube(0.5, 0.02, E, G, rho);
    s.GAy = s.GAz = 0.0;
    s.rIy = s.rIz = 0.0;
    model m;
    m.gravity = Vec3::Zero();
    const int tip = cantilever(m, L, 20, s);
    point_load &pl = static_cast<point_load&>(m.add_load(std::make_unique<point_load>(tip, Vec3(0, 0, -1000.0), Vec3::Zero())));
    static_solve(m, 0.0);
    pl.F.setZero();
    timeint_settings ts;
    ts.rho_inf = 1.0;
    timeint ti(ts);
    ti.initialize(m, 0.0);
    const double w1 = 1.875104*1.875104*std::sqrt(s.EIy/(s.m*L*L*L*L));
    const double T1 = 2.0*M_PI/w1, h = T1/200.0;
    double zprev = m.bodies[tip].x(2);
    std::vector<double> zc;
    for(int i=0; i<1200; ++i)
    {
        ti.step(m, i*h, h);
        const double z = m.bodies[tip].x(2);
        if(zprev*z < 0.0) zc.push_back(i*h + h*zprev/(zprev - z));
        zprev = z;
    }
    const double T = zc.size() > 2 ? 2.0*(zc.back() - zc.front())/(zc.size() - 1) : 0.0;
    check(std::fabs(T - T1) < 0.01*T1, fmt("first period %.5f s vs Euler-Bernoulli %.5f s", T, T1));
}

static void test_rigid_motion()
{
    std::cout<<"free beam in rigid motion"<<std::endl;
    beam_section s = beam_section::tube(0.3, 0.01, E, G, rho);
    model m;
    m.gravity = Vec3::Zero();
    frame_builder fb(m);
    member_spec ms;
    ms.name = "f"; ms.p1 = Vec3(-2, 0, 0); ms.p2 = Vec3(2, 0, 0); ms.n = 4; ms.s1 = s;
    const std::vector<int> st = fb.add_member(ms);
    fb.finish();
    // rigid rotation about z with 0.5 rad/s plus translation
    const Vec3 W(0.0, 0.0, 0.5), V(0.3, 0.1, 0.0);
    for(int i : st)
    {
        rigid_body &b = m.bodies[i];
        b.v = V + W.cross(b.x);
        b.w = b.R.transpose()*W;
    }
    timeint_settings ts;
    ts.rho_inf = 1.0;
    timeint ti(ts);
    ti.initialize(m, 0.0);
    const double E0 = m.kinetic_energy();
    double Umax = 0.0;
    for(int i=0; i<400; ++i)
    {
        ti.step(m, i*0.01, 0.01);
        Umax = std::max(Umax, m.potential_energy());
    }
    check(Umax < 1.0e-6*E0, fmt("strain energy stays below %.1e of the kinetic energy", Umax/E0));
    check(std::fabs(m.kinetic_energy() + m.potential_energy() - E0) < 1.0e-4*E0, fmt("energy kept (rel. error %.1e)", std::fabs(m.kinetic_energy() + m.potential_energy() - E0)/E0));
}

static void test_body_end()
{
    std::cout<<"beam end on a rigid body at an offset"<<std::endl;
    // rigid block at the origin held by a fixed joint to the ground; a cantilever starts 1 m off its CoG
    const double L = 5.0;
    const beam_section s = beam_section::tube(0.3, 0.01, E, G, rho);
    model m;
    m.gravity = Vec3::Zero();
    rigid_body blk;
    blk.name = "block"; blk.mass = 100.0; blk.inertia = Vec3(10, 10, 10).asDiagonal();
    m.add_body(blk);
    m.add_constraint(std::make_unique<joint>(m, joint_type::fixed, -1, 0, Vec3::Zero(), Vec3::UnitZ()));
    frame_builder fb(m);
    member_spec ms;
    ms.name = "arm"; ms.p1 = Vec3(1, 0, 0); ms.p2 = Vec3(1 + L, 0, 0); ms.n = 5; ms.s1 = s; ms.body1 = 0;
    const std::vector<int> st = fb.add_member(ms);
    fb.finish();
    const double P = 50.0;
    m.add_load(std::make_unique<point_load>(st.back(), Vec3(0, 0, -P), Vec3::Zero()));
    timeint ti;
    ti.initialize(m, 0.0);
    // static: the joint carries P and the moment P (L + 1)
    const static_result r = static_solve(m, 0.0);
    const double w = -m.bodies[st.back()].x(2);
    check(r.converged && std::fabs(w - P*L*L*L/(3.0*s.EIy) - P*L/s.GAz) < 1.0e-5*w, fmt("deflection %.6e m as a cantilever clamped at the offset", w));
}

int main()
{
    test_cantilever();
    test_rollup();
    test_bend45();
    test_frequency();
    test_rigid_motion();
    test_body_end();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
