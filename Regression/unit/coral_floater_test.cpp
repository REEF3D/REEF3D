// REEF3D::CORAL unit test of strip-theory hydrodynamics with nonlinear hydrostatics (coral/hydro/coral_strip_hydro): the
// buoyancy of a vertical column at any draft and its heave stiffness rho g A; the pitch stiffness of two columns
// rho g (sum A r^2 + V z_B) of the tilted body; the Morison inertia force on a held column in deep-water waves against the
// integral of the linear kinematics; the dynamic pressure on the end face against linear theory; the free heave of a
// column with its end added mass against T = 2 pi sqrt((m + Ca rho D^3/6)/(rho g A)); kinematics cached at the
// undisplaced strips = exact on a held body; the keyword.
// Architect: Hans Bihs
#include"coral/hydro/coral_strip_hydro.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_gravity.h"
#include"coral/core/coral_timeint.h"
#include"coral/io/coral_input.h"
#include<cmath>
#include<cstdio>
#include<iostream>
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
static const double RHO = 1025.0, G = 9.81;

static int body(model &m, const Vec3 &x, double mass, const Vec3 &I)
{
    rigid_body b;
    b.mass = mass; b.inertia = I.asDiagonal(); b.x = x;
    return m.add_body(b);
}

static void test_hydrostatics()
{
    std::cout<<"hydrostatics"<<std::endl;
    model m;
    const int b = body(m, Vec3(0, 0, -5.0), 1e6, Vec3(1e7, 1e7, 1e7));
    strip_member c;
    c.p1 = Vec3(0, 0, -10.0); c.p2 = Vec3(0, 0, 5.0); c.D = 6.0; c.n = 30; c.end1 = true;
    strip_hydro s(m, b, {c}, nullptr, RHO, G, Vec3::Zero(), "col");
    const double A = M_PI*9.0;
    Vec3 F, M; double V;
    double e = 0.0;
    for(double dz : {0.0, 0.3, -0.77, 1.5})
    {
        m.bodies[b].x = Vec3(0, 0, -5.0 + dz);
        s.total(m, 0.0, F, M, V);
        e = std::max(e, std::fabs(F(2) - RHO*G*A*(10.0 - dz))/(RHO*G*A*10.0));
    }
    check(e < 1e-12, fmt("buoyancy rho g A (draft) of a column at any draft, between strips too (err %.1e)", e));
    m.bodies[b].x = Vec3(0, 0, -5.0);
    s.total(m, 0.0, F, M, V);
    const double F0 = F(2);
    m.bodies[b].x(2) += 0.01;
    s.total(m, 0.0, F, M, V);
    check(std::fabs((F0 - F(2))/0.01/(RHO*G*A) - 1.0) < 1e-9, "heave stiffness rho g A");
    // two columns at +-r: pitch stiffness about the CoG at z_G of the hydrostatic loads
    model m2;
    const double r = 20.0, zG = -6.0;
    const int b2 = body(m2, Vec3(0, 0, zG), 1e6, Vec3(1e8, 1e8, 1e8));
    strip_member c1 = c, c2 = c;
    c1.p1 = Vec3(r, 0, -10.0); c1.p2 = Vec3(r, 0, 5.0);
    c2.p1 = Vec3(-r, 0, -10.0); c2.p2 = Vec3(-r, 0, 5.0);
    strip_hydro s2(m2, b2, {c1, c2}, nullptr, RHO, G, Vec3::Zero(), "cols");
    const double th = 1.0e-4;
    m2.bodies[b2].R = exp_so3(Vec3(0, th, 0));
    s2.total(m2, 0.0, F, M, V);
    const double Vt = 2*A*10.0, zB = -5.0;
    const double C55 = RHO*G*(2*A*r*r + Vt*(zB - zG));
    check(std::fabs(-M(1)/th/C55 - 1.0) < 2e-3, fmt("pitch stiffness of two columns at 20 m = rho g (sum A r^2 + V (z_B - z_G)) = %.4g N m/rad (%.2f %%)",
                                                     C55, (-M(1)/th/C55 - 1.0)*100));
}

static void test_waves()
{
    std::cout<<"waves"<<std::endl;
    auto w = std::make_shared<wave_field>(-1.0, G);
    const double a = 1.0, T = 10.0, om = 2*M_PI/T, k = om*om/G;
    w->add_component(a, om);
    model m;
    const int b = body(m, Vec3(0, 0, -20.0), 1e6, Vec3(1e7, 1e7, 1e7));
    m.add_constraint(std::make_unique<joint>(m, joint_type::fixed, -1, b, Vec3(0, 0, -20.0), Vec3::UnitZ()));
    strip_member c;
    c.p1 = Vec3(0, 0, -40.0); c.p2 = Vec3(0, 0, 5.0); c.D = 2.0; c.n = 400; c.end1 = true; c.Ca = 1.0; c.Cm = 2.0; c.Cd = 0.0;
    strip_hydro s(m, b, {c}, w, RHO, G, Vec3::Zero(), "col");
    const double A = M_PI;
    // t = T/4: crest passed, u = 0 at the surface in phase with the acceleration's maximum: eta = a cos(om t) = 0
    double e = 0.0, emax = 0.0;
    for(double t : {0.0, 2.5, 4.0})
    {
        Vec3 F, M; double V;
        s.total(m, t, F, M, V);
        // the integral of rho A Cm a_x over the wet length, linear kinematics (z = 0 above the still water level)
        const double eta = a*std::cos(om*t);
        double Fx = 0.0;
        const int N = 20000;
        for(int i=0; i<N; ++i)
        {
            const double z = -40.0 + (i + 0.5)/N*45.0;
            if(z > eta) continue;
            const double zz = std::min(z, 0.0);
            Fx -= RHO*A*2.0*om*om*a*std::exp(k*zz)*std::sin(om*t)*45.0/N;
        }
        e = std::max(e, std::fabs(F(0) - Fx));
        emax = std::max(emax, std::fabs(Fx));
        // dynamic pressure at the end face (-40 m): rho g a e^(k z) cos(om t), upwards on the bottom, and the
        // added mass of the end, Ca_end rho D^3/6 times the vertical acceleration -om^2 a e^(k z) cos(om t)
        const double pz = (RHO*G*A - RHO*8.0/6.0*om*om)*a*std::exp(-k*40.0)*std::cos(om*t);
        // added to the buoyancy and the normal Morison force (no vertical component): F_z - buoyancy
        const double buoy = RHO*G*A*(std::min(eta, 5.0) + 40.0);
        e = std::max(e, std::fabs(F(2) - buoy - pz)*1e-3);
    }
    check(e < 1e-3*emax, fmt("held column in deep-water waves: Morison inertia force = the integral of the linear kinematics, "
                             "end face: pressure rho g a e^(kz) and added mass (largest difference %.2e of %.3e N)", e, emax));
    // cached kinematics at the undisplaced strips = exact on a held body
    strip_hydro x(m, b, {c}, w, RHO, G, Vec3::Zero(), "colx");
    x.exact_kinematics = true;
    Vec3 F1, M1, F2, M2; double V1, V2;
    s.total(m, 3.3, F1, M1, V1); x.total(m, 3.3, F2, M2, V2);
    check((F1 - F2).norm() < 1e-9*F1.norm() && (M1 - M2).norm() < 1e-9*M1.norm(), "kinematics cached at the undisplaced strips = exact on a held body");
}

static void test_heave_decay()
{
    std::cout<<"free heave"<<std::endl;
    // a column D 6 m, draft 10 m, floating freely in heave (prismatic joint), mass rho V
    model m;
    const double D = 6.0, A = M_PI*D*D/4, V = A*10.0, mass = RHO*V;
    const int b = body(m, Vec3(0, 0, -5.0), mass, Vec3(1e7, 1e7, 1e7));
    m.gravity = Vec3(0, 0, -G);
    m.add_load(std::make_unique<gravity_load>());
    m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, b, Vec3(0, 0, -5.0), Vec3::UnitZ()));
    strip_member c;
    c.p1 = Vec3(0, 0, -10.0); c.p2 = Vec3(0, 0, 5.0); c.D = D; c.n = 30; c.end1 = true; c.Ca_end = 1.0; c.Cd_end = 0.0; c.Cd = 0.0;
    m.add_load(std::make_unique<strip_hydro>(m, b, std::vector<strip_member>{c}, nullptr, RHO, G, Vec3::Zero(), "col"));
    m.bodies[b].x(2) += 0.2;
    timeint_settings ts; ts.rho_inf = 1.0;
    timeint ti(ts);
    ti.initialize(m, 0.0);
    const double dt = 0.01;
    std::vector<double> z;
    for(int i=0; i<4000; ++i) {ti.step(m, i*dt, dt); z.push_back(m.body(b).x(2) + 5.0);}
    double t0 = -1, t1 = -1; int n = 0;
    for(size_t i=1; i<z.size(); ++i)
        if(z[i-1] < 0 && z[i] >= 0)
        {
            const double tz = (i - 1 + z[i-1]/(z[i-1] - z[i]))*dt + dt;
            if(t0 < 0) t0 = tz;
            t1 = tz; ++n;
        }
    const double Tm = (t1 - t0)/(n - 1);
    const double Te = 2*M_PI*std::sqrt((mass + RHO*D*D*D/6.0)/(RHO*G*A));
    check(std::fabs(Tm/Te - 1.0) < 2e-3, fmt("heave period %.4f s = 2 pi sqrt((m + rho D^3/6)/(rho g A)) = %.4f s", Tm, Te));
}

static void test_input()
{
    std::cout<<"keyword"<<std::endl;
    coral_case c;
    c.read_string("time dt 0.05 end 1\nbody f mass 1e6 cog 0 0 -5\n"
                  "strip_member col body f p1 0 0 -10 p2 0 0 5 D 6 Cd 0.8 n 15 end1 Ca_end 0.8 exact_kinematics\n", ".");
    load *l = nullptr;
    for(auto &q : c.m.load_list) if(q->name()=="col") l = q.get();
    std::vector<std::string> n; std::vector<double> v;
    if(l) {l->channel_names(n); l->channel_values(c.m, 0.0, v);}
    check(l && n.size()==7 && n[6]=="col_V" && std::fabs(v[6] - M_PI*9.0*10.0) < 1e-9 && std::fabs(v[2] - 1025*9.81*M_PI*9.0*10.0) < 1e-6,
          "strip_member: channels _Fx .. _Mz, _V; submerged volume and buoyancy");
}

int main()
{
    test_hydrostatics();
    test_waves();
    test_heave_decay();
    test_input();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
