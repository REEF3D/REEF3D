// Standalone tests of the REEF3D::CORAL Craig-Bampton superelement (src/coral/structure) and of the
// generalized coordinates of the core: no MPI, no REEF3D.
// Build:  see coral_tests.sh (links all CORAL objects)
// Run:    ./coral_superelement_test
//
// Tests: element matrices (rigid-body modes free of strain energy, mass of the element); cantilever
// clamped at the ground with a tip body: static tip deflection (Timoshenko) and deflection under
// self-weight (q L^4/8EI) with all modes kept; fixed-interface frequencies against the clamped-clamped
// beam; the first two frequencies of the cantilever with a tip mass against the Euler-Bernoulli
// frequency equation with few modes; free vibration in the time domain (period, Newton iterations); a jacket-like
// tripod reduced to a superelement (10 modes) under a frame tower against the same structure as a
// full frame (four times finer, lumped mass).
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_gravity.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_static.h"
#include"coral/structure/coral_frame.h"
#include"coral/structure/coral_superelement.h"
#include"coral/linear/coral_modes.h"
#include<cmath>
#include<cstdio>
#include<iostream>
#include<memory>
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
    char s[256];
    std::snprintf(s, sizeof(s), f, a, b, c, d);
    return s;
}

class point_load : public load
{
public:
    point_load(int b, const Vec3 &F) : b_(b), F(F) {}
    std::string name() const override {return "point_load";}
    coupling_kind coupling() const override {return coupling_kind::local;}
    void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
    void eval(const model &m, double, VecX &f) const override {m.add_force(f, b_, m.body(b_).x, F);}
    int b_;
    Vec3 F;
};

static const double E = 2.1e11, G = 8.08e10, rho = 7850.0;

static void test_matrices()
{
    std::cout<<"element matrices"<<std::endl;
    beam_section s = beam_section::tube(0.6, 0.02, E, G, rho);
    const double L = 3.0;
    const Eigen::Matrix<double,12,12> K = beam_stiffness_local(s, L), M = beam_mass_local(s, L);
    // rigid-body motions: translations and rotations about the first node
    double emax = 0.0;
    for(int r=0; r<6; ++r)
    {
        Eigen::Matrix<double,12,1> u = Eigen::Matrix<double,12,1>::Zero();
        if(r < 3) {u(r) = 1.0; u(6+r) = 1.0;}
        else
        {
            const Vec3 th = Vec3::Unit(r-3), rb = Vec3(L, 0, 0);
            u.segment<3>(3) = th; u.segment<3>(9) = th;
            u.segment<3>(6) = th.cross(rb);
        }
        emax = std::max(emax, (K*u).norm()/K.norm());
    }
    check(emax < 1.0e-12, fmt("rigid-body motions carry no force (%.1e)", emax));
    Eigen::Matrix<double,12,1> ux = Eigen::Matrix<double,12,1>::Zero(), uy = ux;
    ux(0) = ux(6) = 1.0;
    uy(1) = uy(7) = 1.0;
    const double mL = s.m*L;
    check(std::fabs(ux.dot(M*ux) - mL) < 1.0e-12*mL && std::fabs(uy.dot(M*uy) - mL - 0.0) < 1.0e-12*mL,
          fmt("translation mass m L = %.3f kg in x and y", mL));
}

// cantilever along x clamped at the ground, tip node on a small rigid body; n elements, nmode modes
static superelement &cantilever(model &m, double L, int n, const beam_section &s, int nmode, double Mtip, const Vec3 &g)
{
    rigid_body tip;
    tip.name = "tip";
    tip.mass = Mtip;
    tip.inertia = 1.0e-3*Mat3::Identity();
    tip.x = Vec3(L, 0, 0);
    m.add_body(tip);
    member_spec ms;
    ms.name = "c"; ms.p2 = Vec3(L, 0, 0); ms.n = n; ms.s1 = s; ms.body1 = -1; ms.body2 = 0;
    auto se = std::make_unique<superelement>(m, "se", std::vector<member_spec>{ms}, nmode, 0.0, 0.0, g);
    superelement &r = *se;
    m.add_load(std::move(se));
    return r;
}

static void test_static()
{
    std::cout<<"cantilever: static"<<std::endl;
    const double L = 8.0, P = 1.0e4;
    const beam_section s = beam_section::tube(0.5, 0.02, E, G, rho);
    {
        model m;
        m.gravity = Vec3::Zero();
        cantilever(m, L, 10, s, 5, 1.0, Vec3::Zero());
        m.add_load(std::make_unique<point_load>(0, Vec3(0, 0, -P)));
        const static_result r = static_solve(m, 0.0);
        const double w = -(m.bodies[0].x(2));
        const double we = P*L*L*L/(3.0*s.EIy) + P*L/s.GAz;
        check(r.converged && std::fabs(w - we) < 1.0e-6*we, fmt("tip load: %.6e m vs Timoshenko %.6e m", w, we));
    }
    {
        // self-weight with all modes: exact static condensation
        model m;
        beam_section s2 = s;
        s2.GAy = s2.GAz = 0.0;
        superelement &se = cantilever(m, L, 10, s2, 1000, 1.0e-6, m.gravity);
        // the tip body itself has (almost) no weight; no gravity load needed for it
        const static_result r = static_solve(m, 0.0);
        const double w = -(m.bodies[0].x(2));
        const double we = s2.m*9.81*L*L*L*L/(8.0*s2.EIy);
        check(r.converged && se.nmode()==54 && std::fabs(w - we) < 1.0e-6*we, fmt("self-weight, all %.0f modes: %.6e m vs q L^4/8EI %.6e m", se.nmode(), w, we));
    }
}

static void test_frequencies()
{
    std::cout<<"cantilever: frequencies"<<std::endl;
    const double L = 8.0;
    beam_section s = beam_section::tube(0.5, 0.02, E, G, rho);
    s.GAy = s.GAz = 0.0;
    s.rIy = s.rIz = 0.0;
    const double c = std::sqrt(s.EIy/(s.m*L*L*L*L))/(2.0*M_PI);
    {
        model m;
        superelement &se = cantilever(m, L, 20, s, 4, 1.0, Vec3::Zero());
        // fixed-interface: clamped at both ends (the tip body is the interface)
        const double fcc = 4.730041*4.730041*c;
        check(std::fabs(se.frequencies()(0) - fcc) < 1.0e-4*fcc, fmt("first fixed-interface mode %.5f Hz vs clamped-clamped %.5f Hz", se.frequencies()(0), fcc));
    }
    {
        // tip mass equal to the beam mass: lambda from 1 + cos cosh + mu lambda (cos sinh - sin cosh) = 0
        const double mu = 1.0;
        auto fe = [&](double l) {return 1.0 + std::cos(l)*std::cosh(l) + mu*l*(std::cos(l)*std::sinh(l) - std::sin(l)*std::cosh(l));};
        double lam[2];
        int k = 0;
        for(double l=0.1; l<10.0 && k<2; l+=0.01)
            if(fe(l)*fe(l + 0.01) < 0.0)
            {
                double a = l, b = l + 0.01;
                for(int it=0; it<60; ++it) {const double cc = 0.5*(a + b); if(fe(a)*fe(cc) <= 0.0) b = cc; else a = cc;}
                lam[k++] = 0.5*(a + b);
            }
        model m;
        m.gravity = Vec3::Zero();
        cantilever(m, L, 20, s, 3, mu*s.m*L, Vec3::Zero());
        m.bodies[0].inertia = 1.0e-6*Mat3::Identity();
        const mode_result r = modes(m, 0.0, 12);
        // bending pairs (y and z): take distinct frequencies
        std::vector<double> f;
        for(int q=0; q<r.frequency.size(); ++q)
            if(f.empty() || std::fabs(r.frequency(q) - f.back()) > 1.0e-6*r.frequency(q)) f.push_back(r.frequency(q));
        const double f1 = lam[0]*lam[0]*c, f2 = lam[1]*lam[1]*c;
        auto best = [&](double fx) {double b = 1.0e30; for(double x : f) b = std::min(b, std::fabs(x - fx)/fx); return b;};
        check(best(f1) < 1.0e-4 && best(f2) < 1.0e-3, fmt("tip mass = beam mass, 3 modes: f1 %.1e, f2 %.1e relative to the frequency equation (%.4f, %.4f Hz)", best(f1), best(f2), f1, f2));
    }
}

static void test_dynamic()
{
    std::cout<<"cantilever: free vibration"<<std::endl;
    const double L = 8.0;
    beam_section s = beam_section::tube(0.5, 0.02, E, G, rho);
    model m;
    m.gravity = Vec3::Zero();
    cantilever(m, L, 10, s, 4, 500.0, Vec3::Zero());
    m.bodies[0].inertia = 1.0*Mat3::Identity();
    const mode_result r = modes(m, 0.0, 1);
    // push the tip, release
    auto pl = std::make_unique<point_load>(0, Vec3(0, 0, -2.0e4));
    point_load *p = pl.get();
    m.add_load(std::move(pl));
    static_solve(m, 0.0);
    p->F.setZero();
    timeint_settings ts;
    ts.rho_inf = 1.0;
    timeint ti(ts);
    ti.initialize(m, 0.0);
    const double E0 = m.kinetic_energy() + m.potential_energy();
    const double T = 1.0/r.frequency(0), h = T/100.0;
    double zp = m.bodies[0].x(2), t = 0.0, tc[3];
    int nc = 0;
    double ke_modes = 0.0;
    for(int i=0; i<500 && nc<3; ++i)
    {
        ti.step(m, t, h);
        t += h;
        const double z = m.bodies[0].x(2);
        if(zp < 0.0 && z >= 0.0) tc[nc++] = t - h*z/(z - zp);
        zp = z;
    }
    (void)ke_modes;
    (void)E0;
    // light tip: the superelement's own mass dominates the Newton matrix
    model ml;
    ml.gravity = Vec3::Zero();
    cantilever(ml, L, 10, s, 4, 1.0, Vec3::Zero());
    ml.bodies[0].v = Vec3(0.0, 0.0, 0.5);
    timeint tl;           // rho_inf 0.8: (1 - alpha_f)/(1 - alpha_m) != 1
    tl.initialize(ml, 0.0);
    for(int i=0; i<200; ++i) tl.step(ml, i*h, h);
    const double Tm = nc==3 ? tc[2] - tc[1] : 0.0;
    check(nc==3 && std::fabs(Tm - T) < 2.0e-3*T, fmt("period %.5f s vs first mode %.5f s", Tm, T));
    // with the exact tangent (mass, stiffness) Newton needs about three iterations per step at tol 1e-12
    const double ips = (double)tl.stats.iterations/tl.stats.steps;
    check(ips < 3.5, fmt("%.2f Newton iterations per step", ips));
}

// tripod-like substructure: three legs from the seabed to a transition piece (rigid body), braced;
// a tower of frame members on the transition piece with a top mass
static void tripod(model &m, bool super, int nmode, int nf)
{
    rigid_body tp;
    tp.name = "tp"; tp.mass = 2.0e5; tp.inertia = Vec3(4e6, 4e6, 6e6).asDiagonal(); tp.x = Vec3(0, 0, 20.0);
    m.add_body(tp);
    const beam_section leg = beam_section::tube(1.2, 0.04, E, G, rho), br = beam_section::tube(0.6, 0.02, E, G, rho);
    std::vector<member_spec> sub;
    const double R = 8.0;
    Vec3 base[3], top[3];
    for(int i=0; i<3; ++i)
    {
        const double a = 2.0*M_PI*i/3.0;
        base[i] = Vec3(R*std::cos(a), R*std::sin(a), 0.0);
        top[i] = Vec3(3.0*std::cos(a), 3.0*std::sin(a), 20.0);
    }
    for(int i=0; i<3; ++i)
    {
        member_spec l;
        l.name = "leg" + std::to_string(i); l.p1 = base[i]; l.p2 = top[i]; l.n = 6*nf; l.s1 = leg; l.body1 = -1; l.body2 = 0;
        l.orient = Vec3::UnitX();
        sub.push_back(l);
        member_spec b;
        b.name = "br" + std::to_string(i); b.p1 = base[i] + 0.5*(top[i] - base[i]); b.p2 = base[(i+1)%3] + 0.5*(top[(i+1)%3] - base[(i+1)%3]);
        b.n = 4*nf; b.s1 = br;
        sub.push_back(b);
    }
    // the brace ends meet the legs at their mid-nodes (n = 6: station 3)
    frame_builder fb(m);
    if(super)
        m.add_load(std::make_unique<superelement>(m, "tripod", sub, nmode, 0.0, 0.0, Vec3::Zero()));
    else
        for(const auto &s : sub) fb.add_member(s);
    member_spec tw;
    tw.name = "tower"; tw.p1 = Vec3(0, 0, 20.0); tw.p2 = Vec3(0, 0, 80.0); tw.n = 12; tw.s1 = beam_section::tube(4.0, 0.03, E, G, rho);
    tw.body1 = 0; tw.orient = Vec3::UnitX();
    const std::vector<int> st = fb.add_member(tw);
    fb.finish();
    m.bodies[st.back()].mass += 3.0e5;
}

static void test_tripod()
{
    std::cout<<"tripod superelement under a frame tower against the full frame"<<std::endl;
    model mf, ms;
    mf.gravity = ms.gravity = Vec3::Zero();
    // the frame (lumped mass) needs a finer mesh than the superelement (consistent mass)
    tripod(mf, false, 0, 4);
    tripod(ms, true, 10, 1);
    const mode_result rf = modes(mf, 0.0, 6), rs = modes(ms, 0.0, 6);
    double e = 0.0;
    for(int q=0; q<6; ++q) e = std::max(e, std::fabs(rs.frequency(q) - rf.frequency(q))/rf.frequency(q));
    check(rf.ok && rs.ok && e < 4.0e-3, fmt("first six modes within %.2f%% (f1 %.4f Hz, f3 %.4f Hz); %.0f DOF in the superelement model", 100.0*e, rs.frequency(0), rs.frequency(2), ms.ndof()));
}

int main()
{
    test_matrices();
    test_static();
    test_frequencies();
    test_dynamic();
    test_tripod();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
