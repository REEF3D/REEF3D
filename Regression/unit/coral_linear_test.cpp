// Standalone tests of the REEF3D::CORAL linear analysis (src/coral/linear): no MPI, no REEF3D.
// Build:  see coral_tests.sh (links all CORAL objects)
// Run:    ./coral_linear_test
//
// Tests: spring-mass (sqrt(k/m)), damped oscillator (damping ratio c/(2 sqrt(k m))), pendulum on a
// revolute joint (geometric stiffness of the constraint force: sqrt(m g L/(m L^2 + J))), cantilever
// beam: first three bending frequencies (Euler-Bernoulli), first torsion and axial frequencies,
// mode shape of the first bending mode; a mass on top of a cantilever against the Rayleigh estimate;
// the NREL 5 MW tower (Euler-Bernoulli sections, gravity, static equilibrium first) with a point mass
// and with the rigid rotor-nacelle assembly against a planar finite-element reference (validation 02).
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_gravity.h"
#include"coral/core/coral_constraint.h"
#include"coral/structure/coral_linear6.h"
#include"coral/structure/coral_frame.h"
#include"coral/core/coral_static.h"
#include"coral/linear/coral_modes.h"
#include"coral/io/coral_modes_out.h"
#include<algorithm>
#include<array>
#include<cmath>
#include<vector>
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

static void test_oscillator()
{
    std::cout<<"oscillators"<<std::endl;
    model m;
    m.gravity = Vec3::Zero();
    rigid_body b;
    b.mass = 4.0;
    b.inertia = Vec3(1, 2, 3).asDiagonal();
    m.add_body(b);
    Mat6 K = Mat6::Zero(), C = Mat6::Zero();
    K(0,0) = 400.0; K(1,1) = 100.0; K(2,2) = 900.0; K(3,3) = 10.0; K(4,4) = 20.0; K(5,5) = 30.0;
    C(2,2) = 6.0;
    m.add_load(std::make_unique<linear6>(m, 0, Vec3::Zero(), K, C));
    const mode_result u = modes(m, 0.0, 6, false);
    const double ex[6] = {std::sqrt(10.0/1.0), std::sqrt(100.0/4.0), std::sqrt(20.0/2.0), std::sqrt(30.0/3.0), std::sqrt(400.0/4.0), std::sqrt(900.0/4.0)};
    std::vector<double> e(ex, ex+6);
    std::sort(e.begin(), e.end());
    double err = 0.0;
    for(int q=0; q<6; ++q) err = std::max(err, std::fabs(2.0*M_PI*u.frequency(q) - e[q])/e[q]);
    check(u.ok && err < 1.0e-6, fmt("six uncoupled modes sqrt(k/m), sqrt(k/J) (max rel. error %.1e)", err));
    const mode_result d = modes(m, 0.0, 6, true);
    // heave: w = 15, zeta = 6/(2 sqrt(900*4)) = 0.05
    int ih = -1;
    for(int q=0; q<6; ++q) if(std::fabs(2.0*M_PI*d.frequency(q) - 15.0) < 0.1) ih = q;
    check(ih >= 0 && std::fabs(d.damping(ih) - 0.05) < 1.0e-6, fmt("damped heave mode: zeta %.6f (0.05)", ih >= 0 ? d.damping(ih) : -1.0));
}

static void test_pendulum()
{
    std::cout<<"pendulum (constraint force in the stiffness)"<<std::endl;
    const double mass = 2.0, L = 1.5, J = 0.1, g = 9.81;
    model m;
    rigid_body b;
    b.mass = mass;
    b.inertia = Vec3(J, J, J).asDiagonal();
    b.x = Vec3(0, 0, -L);
    m.add_body(b);
    m.add_load(std::make_unique<gravity_load>());
    m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, -1, 0, Vec3::Zero(), Vec3::UnitY()));
    const mode_result r = modes(m, 0.0, 1);
    const double w = std::sqrt(mass*g*L/(mass*L*L + J));
    check(r.ok && std::fabs(2.0*M_PI*r.frequency(0) - w) < 1.0e-5*w, fmt("w %.6f rad/s vs sqrt(m g L/(m L^2 + J)) %.6f rad/s", 2.0*M_PI*r.frequency(0), w));
}

static void test_cantilever()
{
    std::cout<<"cantilever beam modes"<<std::endl;
    const double L = 10.0, E = 2.1e11, G = 8.08e10, rho = 7850.0;
    beam_section s = beam_section::tube(0.5, 0.02, E, G, rho);
    s.GAy = s.GAz = 0.0;     // Euler-Bernoulli
    s.rIy = s.rIz = 0.0;     // no rotary inertia of the section
    model m;
    m.gravity = Vec3::Zero();
    frame_builder fb(m);
    member_spec ms;
    ms.name = "c"; ms.p2 = Vec3(L, 0, 0); ms.n = 30; ms.s1 = s; ms.body1 = -1;
    const std::vector<int> st = fb.add_member(ms);
    fb.finish();
    const mode_result r = modes(m, 0.0, 12);
    const double c = std::sqrt(s.EIy/(s.m*L*L*L*L))/(2.0*M_PI);
    const double bl[3] = {1.875104, 4.694091, 7.854757};
    // bending modes come in pairs (y and z); pick the distinct frequencies
    std::vector<double> fb2;
    for(int q=0; q<r.frequency.size(); ++q)
        if(fb2.empty() || std::fabs(r.frequency(q) - fb2.back()) > 1.0e-6*r.frequency(q))
            fb2.push_back(r.frequency(q));
    for(int k=0; k<3; ++k)
    {
        const double fe = bl[k]*bl[k]*c;
        double best = 1.0e30;
        for(double f : fb2) best = std::min(best, std::fabs(f - fe)/fe);
        check(best < (k==0 ? 3.0e-3 : 1.0e-2), fmt("bending mode %.0f: %.4e relative error to Euler-Bernoulli %.4f Hz", k+1, best, fe));
    }
    const double ft = 0.25*std::sqrt(s.GJ/(s.rIx*L*L));
    const double fa = 0.25*std::sqrt(s.EA/(s.m*L*L));
    double bt = 1.0e30, ba = 1.0e30;
    const mode_result r2 = modes(m, 0.0, 40);
    for(int q=0; q<r2.frequency.size(); ++q)
    {
        bt = std::min(bt, std::fabs(r2.frequency(q) - ft)/ft);
        ba = std::min(ba, std::fabs(r2.frequency(q) - fa)/fa);
    }
    check(bt < 3.0e-3, fmt("first torsion mode: %.1e relative error to %.2f Hz", bt, ft));
    check(ba < 3.0e-3, fmt("first axial mode: %.1e relative error to %.2f Hz", ba, fa));
    // first mode shape: tip transverse deflection largest, root zero
    const VecX &v = r.shape.col(0);
    const int tip = st.back();
    const double tipdef = std::max(std::fabs(v(6*tip+1)), std::fabs(v(6*tip+2)));
    check(std::fabs(tipdef - 1.0) < 1.0e-6, "first mode: unit transverse deflection at the tip");
}

static void test_top_mass()
{
    std::cout<<"cantilever with a top mass"<<std::endl;
    // vertical tower with a heavy point mass at the top: Rayleigh estimate with the static shape
    // w = sqrt(3 EI / (L^3 (M + 33/140 m L)))
    const double L = 20.0, E = 2.1e11, G = 8.08e10, rho = 7850.0, Mtop = 5000.0;
    beam_section s = beam_section::tube(0.8, 0.015, E, G, rho);
    s.GAy = s.GAz = 0.0;
    s.rIy = s.rIz = 0.0;
    model m;
    m.gravity = Vec3::Zero();
    frame_builder fb(m);
    member_spec ms;
    ms.name = "t"; ms.p2 = Vec3(0, 0, L); ms.n = 20; ms.s1 = s; ms.body1 = -1; ms.orient = Vec3::UnitX();
    const std::vector<int> st = fb.add_member(ms);
    fb.finish();
    m.bodies[st.back()].mass += Mtop;
    const mode_result r = modes(m, 0.0, 2);
    const double w = std::sqrt(3.0*s.EIy/(L*L*L*(Mtop + 33.0/140.0*s.m*L)));
    check(std::fabs(2.0*M_PI*r.frequency(0) - w) < 0.01*w, fmt("first mode %.5f Hz vs Rayleigh %.5f Hz", r.frequency(0), w/(2.0*M_PI)));
}

static void test_nrel_tower()
{
    std::cout<<"NREL 5 MW tower with the RNA (Euler-Bernoulli, gravity)"<<std::endl;
    // tower 87.6 m, D 6 -> 3.87 m, wall 0.0351 -> 0.0247 m, 20 elements; RNA 350 t, CoM (-0.2746, 0, 1.9543)
    // above the top. Reference: Dropbox/Claude/CORAL/validation/02_nrel5mw_tower/reference.py (200 elements)
    const double H = 87.6, E = 210e9, G = 80.8e9, rho = 8500.0, M = 350000.0;
    const double xc = -0.27457428571428570, zc = 1.9542857142857142;
    for(int variant=0; variant<2; ++variant)
    {
        model m;
        frame_builder fb(m);
        std::vector<int> top;
        for(int i=0; i<20; ++i)
        {
            const double f = (i + 0.5)/20.0, D = 6.0 - 2.13*f, t = 0.0351 - 0.0104*f;
            beam_section s = beam_section::tube(D, t, E, G, rho);
            s.GAy = s.GAz = 0.0;
            s.rIy = s.rIz = 0.0;
            member_spec ms;
            ms.name = "t" + std::to_string(i+1);
            ms.p1 = Vec3(0, 0, H*i/20.0); ms.p2 = Vec3(0, 0, H*(i+1)/20.0); ms.n = 1; ms.s1 = s;
            ms.orient = Vec3::UnitX(); ms.body1 = i==0 ? -1 : -2;
            top = fb.add_member(ms);
        }
        fb.finish();
        if(variant==0)
            m.bodies[top.back()].mass += M;
        else
        {
            rigid_body rna;
            rna.mass = M;
            rna.inertia << 41398994.571428575, 0, 339233.58857142855, 0, 25630439.837097142, 0, 339233.58857142855, 0, 25598571.26566857;
            rna.x = Vec3(xc, 0, H + zc);
            const int b = m.add_body(rna);
            m.add_constraint(std::make_unique<joint>(m, joint_type::fixed, top.back(), b, Vec3(0, 0, H), Vec3::UnitZ()));
        }
        m.add_load(std::make_unique<gravity_load>());
        const static_result st = static_solve(m, 0.0);
        const mode_result r = modes(m, 0.0, 8);
        std::vector<double> fa, ss;
        for(int q=0; q<r.frequency.size(); ++q)
        {
            const std::array<double,6> e = energy_split(m, r, q);
            if(e[0] + e[4] > 0.5) fa.push_back(r.frequency(q));
            else if(e[1] + e[3] > 0.5) ss.push_back(r.frequency(q));
        }
        const double ref[2][4] = {{0.33048122, 0.33048122, 3.06418737, 3.06418737}, {0.31562006, 0.31329227, 2.18777222, 1.90506935}};
        bool ok = st.converged && fa.size() >= 2 && ss.size() >= 2;
        double err[4] = {1, 1, 1, 1};
        if(ok)
        {
            err[0] = fa[0]/ref[variant][0] - 1.0; err[1] = ss[0]/ref[variant][1] - 1.0;
            err[2] = fa[1]/ref[variant][2] - 1.0; err[3] = ss[1]/ref[variant][3] - 1.0;
        }
        for(int k=0; k<4; ++k) ok = ok && std::fabs(err[k]) < (k < 2 ? 1.5e-3 : 3.0e-3);    // 20 lumped-mass elements
        check(ok, fmt(variant==0 ? "point mass: FA1 SS1 FA2 SS2 %+.2e %+.2e %+.2e %+.2e relative to the reference"
                                 : "rigid RNA:  FA1 SS1 FA2 SS2 %+.2e %+.2e %+.2e %+.2e relative to the reference", err[0], err[1], err[2], err[3]));
    }
}

int main()
{
    test_oscillator();
    test_pendulum();
    test_cantilever();
    test_top_mass();
    test_nrel_tower();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
