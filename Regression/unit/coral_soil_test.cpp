// Standalone tests of the REEF3D::CORAL soil springs (src/coral/soil): no MPI, no REEF3D.
// Build:  see coral_tests.sh (links all CORAL objects)
// Run:    ./coral_soil_test
//
// Tests: API soft clay p-y (pu, 0.5 pu at yc, pu at 8 yc, cyclic 0.72 pu), API sand p-y (wedge and flow
// resistance against the API coefficients C1, C2, C3 for phi = 30 and 35 deg, initial slope k X,
// A pu at large y), t-z and Q-z (Q-z only in compression); spring tangent against finite differences;
// long elastic pile in a linear Winkler soil under a lateral head load and a head moment against
// Hetenyi's semi-infinite beam (y0 = 2 H beta/k, theta0 = 2 H beta^2/k, y0 = 2 M beta^2/k) and under
// an axial load (head stiffness sqrt(EA kv)); a pile in API clay under a large lateral load: the soil
// forces balance the load.
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_static.h"
#include"coral/core/coral_rotation.h"
#include"coral/structure/coral_frame.h"
#include"coral/soil/coral_soil.h"
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
    point_load(int b, const Vec3 &F, const Vec3 &M) : b_(b), F(F), M(M) {}
    std::string name() const override {return "point_load";}
    coupling_kind coupling() const override {return coupling_kind::local;}
    void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
    void eval(const model &m, double, VecX &f) const override {m.add_force(f, b_, m.body(b_).x, F); m.add_moment(f, b_, M);}
    int b_;
    Vec3 F, M;
};

// rotation of a node about global y (small)
static double log_rotation_y(const rigid_body &b) {return -log_so3(b.R)(1);}

static void test_curves()
{
    std::cout<<"API curves"<<std::endl;
    soil_profile clay;
    clay.mudline = 0.0;
    soil_layer c;
    c.type = soil_layer::kind::clay;
    c.top = 0.0; c.bottom = 40.0; c.gamma = 8.0e3; c.su_top = 20.0e3; c.su_bottom = 60.0e3; c.eps50 = 0.01; c.J = 0.5;
    clay.layers.push_back(c);
    const double D = 2.0, X = 5.0;
    const double su = 20.0e3 + 40.0e3*X/40.0, p0 = 8.0e3*X;
    const double pu = std::min(3.0*su + p0 + 0.5*su*X/D, 9.0*su)*D, yc = 2.5*0.01*D;
    const soil_curve py = clay.py(X, D);
    check(std::fabs(py.force(yc) - 0.5*pu) < 1.0e-9*pu && std::fabs(py.force(8.0*yc) - pu) < 1.0e-9*pu && std::fabs(py.force(-20.0*yc) + pu) < 1.0e-9*pu,
          fmt("soft clay static: pu %.4e N/m, 0.5 pu at yc, pu at 8 yc, odd (%.4e N/m at -20 yc)", pu, py.force(-20.0*yc)));
    clay.layers[0].cyclic = true;
    const double XR = 6.0*D/(8.0e3*D/su + 0.5);
    const soil_curve pyc = clay.py(X, D);
    check(X < XR && std::fabs(pyc.force(15.0*yc) - 0.72*pu*X/XR) < 1.0e-9*pu && std::fabs(pyc.force(3.0*yc) - 0.5*std::cbrt(3.0)*pu) < 1.0e-9*pu,
          fmt("soft clay cyclic above XR = %.2f m: 0.72 pu X/XR at 15 yc", XR));

    soil_profile sand;
    soil_layer s;
    s.type = soil_layer::kind::sand;
    s.top = 0.0; s.bottom = 50.0; s.gamma = 10.0e3; s.phi = 30.0; s.k = 16.3e6;
    sand.layers.push_back(s);
    // API coefficients (RP 2A Fig. 6.8.6-1): C1 ~ 1.9, C2 ~ 2.7, C3 ~ 29 at 30 deg; 3.0, 3.4, 54 at 35 deg
    auto pu_sand = [&](double phi, double XX)
    {
        sand.layers[0].phi = phi;
        const soil_curve q = sand.py(XX, D);
        const double A = std::max(0.9, 3.0 - 0.8*XX/D);
        return q.force(1.0e3)/A;
    };
    const double g = 10.0e3;
    const double C3_30 = pu_sand(30.0, 40.0)/(D*g*40.0);           // deep: flow
    const double sh1 = pu_sand(30.0, 1.0), sh2 = pu_sand(30.0, 2.0); // shallow: (C1 X + C2 D) g X
    const double C1_30 = (sh2/(g*2.0) - sh1/(g*1.0))/1.0, C2_30 = (sh1/(g*1.0) - C1_30*1.0)/D;
    const double C3_35 = pu_sand(35.0, 40.0)/(D*g*40.0);
    check(std::fabs(C1_30 - 1.9) < 0.05 && std::fabs(C2_30 - 2.7) < 0.05 && std::fabs(C3_30 - 29.0) < 0.5 && std::fabs(C3_35 - 54.0) < 1.0,
          fmt("sand: C1 %.2f, C2 %.2f, C3 %.1f at 30 deg, C3 %.1f at 35 deg", C1_30, C2_30, C3_30, C3_35));
    sand.layers[0].phi = 30.0;
    const soil_curve ps = sand.py(10.0, D);
    check(std::fabs(ps.slope(0.0)*(ps.y[1]) - ps.p[1]) < 1.0e-12*ps.p[1] && std::fabs(ps.p[1]/ps.y[1] - s.k*10.0) < 0.01*s.k*10.0,
          fmt("sand: initial slope %.4e N/m^2 vs k X %.4e N/m^2 (within the tanh secant)", ps.p[1]/ps.y[1], s.k*10.0));
    // t-z and Q-z
    const soil_curve tz = clay.tz(X, D), qz = clay.qz(X, D);
    const double psi = su/p0, alpha = std::min(1.0, psi <= 1.0 ? 0.5/std::sqrt(psi) : 0.5*std::pow(psi, -0.25));
    check(std::fabs(tz.force(0.01*D) - alpha*su*M_PI*D) < 1.0e-9*su && std::fabs(tz.force(1.0) - 0.9*alpha*su*M_PI*D) < 1.0e-9*su,
          fmt("clay t-z: t_max = alpha su pi D (alpha %.3f), residual 0.9", alpha));
    check(std::fabs(qz.force(0.1*D) - 9.0*su*M_PI*D*D/4.0) < 1.0e-9*su && qz.force(-0.1) == 0.0, "Q-z: 9 su A at 0.1 D, nothing in tension");
}

static void test_tangent()
{
    std::cout<<"spring tangent"<<std::endl;
    model m;
    rigid_body b;
    b.mass = 1.0;
    b.inertia = Mat3::Identity();
    m.add_body(b);
    soil_curve lat, ax, tip;
    lat.y = {0.0, 0.01, 0.05, 0.2}; lat.p = {0.0, 1.0e5, 3.0e5, 4.0e5};
    ax.y = {0.0, 0.005, 0.02}; ax.p = {0.0, 2.0e5, 2.5e5};
    tip = soil_curve::linear(1.0e7); tip.compression_only = true;
    m.add_load(std::make_unique<soil_spring>(m, 0, Vec3(0.1, 0.0, -1.0), lat, ax, tip, "s", 1.5));
    m.bodies[0].x = Vec3(0.013, -0.021, -0.003);
    m.bodies[0].R = exp_so3(Vec3(0.0005, -0.0002, -0.004));
    std::vector<triplet> T;
    m.load_list[0]->tangent(m, 0.0, 1.0, 0.0, 0.0, T);
    Eigen::Matrix<double,6,6> K = Eigen::Matrix<double,6,6>::Zero(), Kf;
    for(const auto &t : T) K(t.row(), t.col()) += t.value();
    const Vec3 x0 = m.bodies[0].x;
    const Mat3 R0 = m.bodies[0].R;
    for(int j=0; j<6; ++j)
    {
        const double e = 1.0e-8;
        VecX fp = VecX::Zero(6), fm = VecX::Zero(6);
        for(int sg=-1; sg<=1; sg+=2)
        {
            m.bodies[0].x = x0; m.bodies[0].R = R0;
            if(j < 3) m.bodies[0].x(j) += sg*e; else m.bodies[0].R = R0*exp_so3(sg*e*Vec3::Unit(j-3));
            m.load_list[0]->eval(m, 0.0, sg < 0 ? fm : fp);
        }
        Kf.col(j) = (fp - fm)/(2.0*e);
    }
    m.bodies[0].x = x0; m.bodies[0].R = R0;
    const double err = (K - Kf).cwiseAbs().maxCoeff()/Kf.cwiseAbs().maxCoeff();
    check(err < 1.0e-6, fmt("lateral, axial, tip and shaft torsion: max difference %.1e relative", err));
}

// vertical tube pile from the mudline (z = 0) down to -L, n elements; returns the stations (head first)
static std::vector<int> pile(model &m, double L, int n, const beam_section &s)
{
    frame_builder fb(m);
    member_spec ms;
    ms.name = "pile"; ms.p1 = Vec3(0, 0, 0); ms.p2 = Vec3(0, 0, -L); ms.n = n; ms.s1 = s; ms.orient = Vec3::UnitX();
    const std::vector<int> st = fb.add_member(ms);
    fb.finish();
    return st;
}

static void test_hetenyi()
{
    std::cout<<"elastic pile in a Winkler soil (Hetenyi)"<<std::endl;
    beam_section s = beam_section::tube(2.0, 0.05, 2.1e11, 8.08e10, 7850.0);
    s.GAy = s.GAz = 0.0;
    const double kh = 5.0e7, kv = 2.0e7, L = 60.0;
    const double beta = std::pow(kh/(4.0*s.EIy), 0.25);
    soil_profile sp;
    sp.mudline = 0.0;
    soil_layer l;
    l.type = soil_layer::kind::linear;
    l.top = 0.0; l.bottom = 100.0; l.kh = kh; l.kv = kv;
    sp.layers.push_back(l);
    for(int c=0; c<3; ++c)
    {
        model m;
        m.gravity = Vec3::Zero();
        const std::vector<int> st = pile(m, L, 240, s);
        add_pile_springs(m, sp, st, 2.0, true, true, false, false, "pile");
        const double H = 1.0e6, M0 = 5.0e6, P = 2.0e7;
        m.add_load(std::make_unique<point_load>(st[0], c==0 ? Vec3(H, 0, 0) : c==2 ? Vec3(0, 0, -P) : Vec3::Zero(),
                                                c==1 ? Vec3(0, M0, 0) : Vec3::Zero()));
        const static_result r = static_solve(m, 0.0);
        const rigid_body &h = m.bodies[st[0]];
        if(c==0)
        {
            const double y0 = 2.0*H*beta/kh, t0 = 2.0*H*beta*beta/kh;
            const double t = -log_rotation_y(h);
            check(r.converged && std::fabs(h.x(0) - y0) < 2.0e-3*y0 && std::fabs(t - t0) < 5.0e-3*t0,
                  fmt("head load: y0 %.5e m vs %.5e m, rotation %.4e vs %.4e", h.x(0), y0, t, t0));
        }
        else if(c==1)
        {
            const double y0 = 2.0*M0*beta*beta/kh;
            check(r.converged && std::fabs(std::fabs(h.x(0)) - y0) < 5.0e-3*y0, fmt("head moment: |y0| %.5e m vs 2 M beta^2/k %.5e m", std::fabs(h.x(0)), y0));
        }
        else
        {
            const double kax = std::sqrt(s.EA*kv)*std::tanh(std::sqrt(kv/s.EA)*L), w = -h.x(2);
            check(r.converged && std::fabs(w - P/kax) < 2.0e-3*P/kax, fmt("axial: settlement %.5e m vs P/(sqrt(EA kv) tanh) %.5e m", w, P/kax));
        }
    }
}

static void test_clay_pile()
{
    std::cout<<"pile in API clay under a large lateral load"<<std::endl;
    beam_section s = beam_section::tube(2.0, 0.05, 2.1e11, 8.08e10, 7850.0);
    soil_profile sp;
    soil_layer c;
    c.type = soil_layer::kind::clay;
    c.top = 0.0; c.bottom = 60.0; c.gamma = 8.0e3; c.su_top = 10.0e3; c.su_bottom = 100.0e3; c.eps50 = 0.01; c.J = 0.5;
    sp.layers.push_back(c);
    model m;
    m.gravity = Vec3::Zero();
    const std::vector<int> st = pile(m, 40.0, 40, s);
    const int ns = add_pile_springs(m, sp, st, 2.0, true, true, true, false, "pile");
    const double H = 3.0e6;
    m.add_load(std::make_unique<point_load>(st[0], Vec3(H, 0, 0), Vec3::Zero()));
    // load steps
    bool ok = true;
    auto *pl = dynamic_cast<point_load*>(m.load_list.back().get());
    for(int k=1; k<=5 && ok; ++k)
    {
        pl->F = Vec3(H*k/5.0, 0, 0);
        ok = static_solve(m, 0.0).converged;
    }
    Vec3 Fs = Vec3::Zero();
    for(const auto &l : m.load_list)
        if(const auto *sp2 = dynamic_cast<const soil_spring*>(l.get())) Fs += sp2->force(m);
    check(ok && ns==40 && std::fabs(Fs(0) + H) < 1.0e-6*H, fmt("load steps converged, %.0f springs below the mudline carry %.6e N of %.6e N", ns, -Fs(0), H));
    check(m.bodies[st[0]].x(0) > 0.0, fmt("head deflection %.4f m", m.bodies[st[0]].x(0)));
}

int main()
{
    test_curves();
    test_tangent();
    test_hetenyi();
    test_clay_pile();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
