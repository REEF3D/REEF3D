// Standalone tests of the REEF3D::CORAL loads on flexible offshore structures: Morison loads on frame
// and superelement members (src/coral/hydro/coral_member_morison), through the material points of the
// core (src/coral/core/coral_points): no MPI, no REEF3D.
// Build:  see coral_tests.sh (links all CORAL objects)
// Run:    ./coral_offshore_test
//
// Tests: buoyancy of a partly submerged column; inertia force on a stiff monopile in a linear regular
// wave against Cm rho A w^2 H/(2k) (the integral of the acceleration over the depth); added mass of a
// submerged cantilever (frequency ratio sqrt(m/(m + Ca rho A))); a cantilever in a current as frame
// members and as a superelement (same static deflection); a point mass at an offset (compound
// pendulum).
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_gravity.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_static.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_points.h"
#include"coral/structure/coral_frame.h"
#include"coral/structure/coral_superelement.h"
#include"coral/hydro/coral_member_morison.h"
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

static const double E = 2.1e11, G = 8.08e10, rho_s = 7850.0, rho = 1025.0, g = 9.81;

// vertical column clamped at z0, top at z1, n elements; Morison on it
static member_morison *column(model &m, double z0, double z1, int n, const beam_section &s, const member_morison_coef &c,
                              std::shared_ptr<const wave_field> w, const Vec3 &current, std::vector<int> *st_out=nullptr)
{
    frame_builder fb(m);
    member_spec ms;
    ms.name = "col"; ms.p1 = Vec3(0, 0, z0); ms.p2 = Vec3(0, 0, z1); ms.n = n; ms.s1 = s; ms.body1 = -1; ms.orient = Vec3::UnitX();
    const std::vector<int> st = fb.add_member(ms);
    fb.finish();
    auto bp = std::make_shared<body_points>();
    std::vector<std::pair<int,int>> el;
    for(int i=0; i<=n; ++i)
    {
        bp->add(m, st[i], ms.p1 + (double)i/n*(ms.p2 - ms.p1));
        if(i > 0) el.emplace_back(i-1, i);
    }
    auto mm = std::make_unique<member_morison>(bp, el, c, w, rho, g, current, "morison");
    member_morison *r = mm.get();
    m.add_load(std::move(mm));
    if(st_out) *st_out = st;
    return r;
}

static void test_buoyancy()
{
    std::cout<<"buoyancy"<<std::endl;
    model m;
    m.gravity = Vec3::Zero();
    member_morison_coef c;
    c.D = 6.0;
    member_morison *mm = column(m, -20.0, 10.0, 30, beam_section::tube(6.0, 0.06, E, G, rho_s), c, nullptr, Vec3::Zero());
    Vec3 F, M;
    mm->total(m, 0.0, F, M);
    const double Fb = rho*g*0.25*M_PI*36.0*20.0;
    check(std::fabs(F(2) - Fb) < 1.0e-9*Fb && F.head<2>().norm() < 1.0e-9*Fb, fmt("submerged 20 m of 30 m: %.6e N vs rho g A 20 m %.6e N", F(2), Fb));
}

static void test_inertia()
{
    std::cout<<"monopile in a regular wave (inertia)"<<std::endl;
    const double d = 20.0, H = 0.2, T = 8.0, D = 6.0;    // small H: the load above the still water level is of second order
    auto w = std::make_shared<wave_field>(d, g);
    w->add_regular(H, T);
    model m;
    m.gravity = Vec3::Zero();
    member_morison_coef c;
    c.D = D; c.Cd = 0.0; c.Ca = 1.0; c.n = 4; c.buoyancy = false;
    // very stiff column: the structure does not move
    beam_section s = beam_section::tube(D, 0.5, 1.0e3*E, 1.0e3*G, rho_s);
    member_morison *mm = column(m, -d, 10.0, 60, s, c, w, Vec3::Zero());
    double Fmax = 0.0;
    for(int i=0; i<400; ++i)
    {
        Vec3 F, M;
        mm->total(m, i*T/400.0, F, M);
        Fmax = std::max(Fmax, std::fabs(F(0)));
    }
    const double om = 2.0*M_PI/T, k = wave_field::wavenumber(om, d, g);
    const double FI = 2.0*rho*0.25*M_PI*D*D*om*om*H/(2.0*k);
    check(std::fabs(Fmax - FI) < 3.0e-3*FI, fmt("max force %.5e N vs Cm rho A w^2 H/(2k) %.5e N", Fmax, FI));
}

static void test_added_mass()
{
    std::cout<<"submerged cantilever: added mass"<<std::endl;
    const double D = 1.0;
    const beam_section s = beam_section::tube(D, 0.02, E, G, rho_s);
    double f[2];
    for(int wet=0; wet<2; ++wet)
    {
        model m;
        m.gravity = Vec3::Zero();
        member_morison_coef c;
        c.D = D; c.Cd = 0.0; c.Ca = wet ? 1.0 : 0.0; c.buoyancy = false; c.n = 4;
        column(m, -40.0, -10.0, 30, s, c, nullptr, Vec3::Zero());
        f[wet] = modes(m, 0.0, 1).frequency(0);
    }
    const double ma = rho*0.25*M_PI*D*D, ratio = std::sqrt(s.m/(s.m + ma));
    check(std::fabs(f[1]/f[0] - ratio) < 2.0e-3*ratio, fmt("f_wet/f_dry %.5f vs sqrt(m/(m + Ca rho A)) %.5f (dry %.4f Hz)", f[1]/f[0], ratio, f[0]));
}

static void test_superelement_current()
{
    std::cout<<"cantilever in a current: frame and superelement"<<std::endl;
    const double D = 1.0;
    const beam_section s = beam_section::tube(D, 0.02, E, G, rho_s);
    const Vec3 U(1.5, 0.0, 0.0);
    member_morison_coef c;
    c.D = D; c.Cd = 1.0; c.Ca = 1.0; c.buoyancy = false; c.n = 2;
    // frame
    model mf;
    mf.gravity = Vec3::Zero();
    std::vector<int> st;
    column(mf, -30.0, -5.0, 25, s, c, nullptr, U, &st);
    // superelement: the column with its top on a small body
    model ms;
    ms.gravity = Vec3::Zero();
    rigid_body top;
    top.name = "top"; top.mass = 1.0; top.inertia = Mat3::Identity(); top.x = Vec3(0, 0, -5.0);
    ms.add_body(top);
    member_spec mp;
    mp.name = "col"; mp.p1 = Vec3(0, 0, -30.0); mp.p2 = Vec3(0, 0, -5.0); mp.n = 25; mp.s1 = s; mp.body1 = -1; mp.body2 = 0; mp.orient = Vec3::UnitX();
    auto se = std::make_unique<superelement>(ms, "se", std::vector<member_spec>{mp}, 20, 0.0, 0.0, Vec3::Zero());
    superelement *sp = se.get();
    ms.add_load(std::move(se));
    std::vector<std::pair<int,int>> el;
    const std::vector<int> nd = sp->member_nodes("col");
    for(size_t i=1; i<nd.size(); ++i) el.emplace_back(nd[i-1], nd[i]);
    ms.add_load(std::make_unique<member_morison>(std::shared_ptr<const point_set>(std::shared_ptr<const point_set>(), sp), el, c, nullptr, rho, g, U, "morison"));
    const static_result rf = static_solve(mf, 0.0), rs = static_solve(ms, 0.0);
    const double xf = mf.bodies[st.back()].x(0), xs = ms.bodies[0].x(0);
    // uniform load q = 1/2 rho Cd D U^2 (the tip moves little): q L^4/(8 EI) + q L^2/(2 GA)
    const double q = 0.5*rho*D*U(0)*U(0), L = 25.0, w = q*L*L*L*L/(8.0*s.EIy) + q*L*L/(2.0*s.GAz);
    check(rf.converged && rs.converged && std::fabs(xs - xf) < 2.0e-3*xf && std::fabs(xf - w) < 1.0e-2*w,
          fmt("tip %.5e m (frame), %.5e m (superelement, 20 modes), q L^4/8EI + q L^2/2GA %.5e m", xf, xs, w));
}

static void test_point_mass()
{
    std::cout<<"point mass at an offset"<<std::endl;
    // massless-ish rod body on a revolute joint at the origin with a point mass 2 m below
    model m;
    rigid_body b;
    b.mass = 1.0; b.inertia = 0.1*Mat3::Identity();
    m.add_body(b);
    const double mp = 50.0, r = 2.0;
    m.add_point_mass(0, Vec3(0, 0, -r), mp, Mat3::Zero());
    m.add_load(std::make_unique<gravity_load>());
    m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, -1, 0, Vec3::Zero(), Vec3::UnitY()));
    const mode_result md = modes(m, 0.0, 1);
    const double w = std::sqrt(mp*9.81*r/(0.1 + mp*r*r));
    check(md.ok && std::fabs(2.0*M_PI*md.frequency(0) - w) < 1.0e-5*w, fmt("compound pendulum %.6f rad/s vs %.6f rad/s", 2.0*M_PI*md.frequency(0), w));
}

int main()
{
    test_buoyancy();
    test_inertia();
    test_added_mass();
    test_superelement_current();
    test_point_mass();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
