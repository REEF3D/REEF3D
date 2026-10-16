// Standalone tests of the REEF3D::CORAL power take-off (src/coral/pto): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         coral_pto_test.cpp ../../src/coral/core/*.cpp ../../src/coral/structure/*.cpp
//         ../../src/coral/hydro/*.cpp ../../src/coral/pto/*.cpp -o coral_pto_test
// Run:    ./coral_pto_test
//
// Tests: linear damper (damping ratio of the decay, absorbed energy = loss of mechanical energy),
// Coulomb friction (stopping time and distance), end stop (penetration v0 sqrt(m/k), rebound),
// latching (held at the turning point), declutching (free motion), power limit, rotational PTO in a
// revolute joint; two-body point absorber (float and spar, PTO in the prismatic joint between them,
// synthetic hydrodynamic database) in regular waves against the frequency domain: relative motion
// and mean absorbed power 1/2 c w^2 |Z1 - Z2|^2.
// PTO element library (shared with the REEF3D 6DOF bodies, X 500 - X 508): Jacobians of the elements
// against finite differences, composite with force and power limits and the declutching window,
// zero crossing interpolation, latching and declutching controllers, hydrodynamic table and B_opt,
// prismatic joint of one body (q, qd, J against finite differences); latching with a stiff latch
// spring in the CORAL joint PTO.
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_gravity.h"
#include"coral/structure/coral_linear6.h"
#include"coral/hydro/coral_hydro_load.h"
#include"coral/pto/coral_pto.h"
#include"coral/pto/coral_pto_elements.h"
#include"coral/pto/coral_pto_joint.h"
#include"coral/pto/coral_pto_controller.h"
#include"coral/core/coral_rotation.h"
#include<fstream>
#include<iostream>
#include<string>
#include<cstdio>
#include<cmath>
#include<complex>
#include<memory>

using namespace coral;
using cplx = std::complex<double>;

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

// body on a prismatic joint along x to the ground, optional spring k to the ground
struct slider
{
    model m;
    joint *j = nullptr;
    pto_load *pto = nullptr;
    slider(double mass, double k, const pto_settings &s, bool spring_in_pto=false)
    {
        rigid_body b;
        b.mass = mass;
        b.inertia = Vec3(1, 1, 1).asDiagonal();
        m.add_body(b);
        j = &static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 0, Vec3::Zero(), Vec3::UnitX())));
        if(k > 0.0 && !spring_in_pto)
        {
            Mat6 K = Mat6::Zero();
            K(0,0) = k;
            m.add_load(std::make_unique<linear6>(m, 0, Vec3::Zero(), K, Mat6::Zero()));
        }
        pto = &static_cast<pto_load&>(m.add_load(std::make_unique<pto_load>(*j, false, s)));
    }
};

static void test_damper()
{
    std::cout<<"linear damper"<<std::endl;
    const double mass = 10.0, k = 1000.0, c = 20.0;
    pto_settings s;
    s.c = c;
    slider sl(mass, k, s);
    sl.m.bodies[0].x(0) = 0.1;
    timeint ti;
    ti.initialize(sl.m, 0.0);
    const double E0 = sl.m.kinetic_energy() + sl.m.potential_energy();
    std::vector<double> peaks;
    double x1 = 0.1, x2 = 0.1;
    for(int i=0; i<3000; ++i)
    {
        ti.step(sl.m, i*1.0e-3, 1.0e-3);
        const double x = sl.m.bodies[0].x(0);
        if(x2 < x1 && x1 > x && x1 > 0.0) peaks.push_back(x1);
        x2 = x1; x1 = x;
    }
    const double E = sl.m.kinetic_energy() + sl.m.potential_energy();
    const double zeta = c/(2.0*std::sqrt(k*mass));
    const double wd = std::sqrt(k/mass)*std::sqrt(1.0 - zeta*zeta);
    const double dec = std::log(peaks[0]/peaks[1]);
    check(std::fabs(dec - 2.0*M_PI*zeta/std::sqrt(1.0 - zeta*zeta)) < 2.0e-3, fmt("logarithmic decrement %.5f vs %.5f", dec, 2.0*M_PI*zeta/std::sqrt(1.0 - zeta*zeta)));
    (void)wd;
    check(std::fabs((E0 - E) - sl.pto->energy()) < 1.0e-4*E0, fmt("absorbed energy %.6f J = loss of mechanical energy %.6f J", sl.pto->energy(), E0 - E));
}

static void test_coulomb()
{
    std::cout<<"Coulomb friction"<<std::endl;
    const double mass = 5.0, Fc = 10.0, v0 = 2.0;
    pto_settings s;
    s.coulomb = Fc;
    s.v_eps = 1.0e-4;
    slider sl(mass, 0.0, s);
    sl.m.bodies[0].v(0) = v0;
    timeint ti;
    ti.initialize(sl.m, 0.0);
    double tstop = -1.0;
    for(int i=0; i<2000; ++i)
    {
        ti.step(sl.m, i*1.0e-3, 1.0e-3);
        if(tstop < 0.0 && sl.m.bodies[0].v(0) < 1.0e-3) tstop = (i+1)*1.0e-3;
    }
    const double xs = 0.5*mass*v0*v0/Fc, ts = mass*v0/Fc;
    check(std::fabs(sl.m.bodies[0].x(0) - xs) < 0.002*xs, fmt("stopping distance %.4f m vs m v0^2/(2 Fc) = %.4f m", sl.m.bodies[0].x(0), xs));
    check(std::fabs(tstop - ts) < 0.005*ts, fmt("stopping time %.4f s vs m v0/Fc = %.4f s", tstop, ts));
}

static void test_end_stop()
{
    std::cout<<"end stop"<<std::endl;
    const double mass = 2.0, v0 = 1.5, ks = 2.0e5;
    pto_settings s;
    s.q_max = 0.5;
    s.k_stop = ks;
    slider sl(mass, 0.0, s);
    sl.m.bodies[0].v(0) = v0;
    timeint ti;
    ti.initialize(sl.m, 0.0);
    double qmax = 0.0;
    for(int i=0; i<2000; ++i)
    {
        ti.step(sl.m, i*2.0e-4, 2.0e-4);
        qmax = std::max(qmax, sl.m.bodies[0].x(0));
    }
    const double pen = v0*std::sqrt(mass/ks);
    check(std::fabs(qmax - 0.5 - pen) < 0.01*pen, fmt("penetration %.5f m vs v0 sqrt(m/k) = %.5f m", qmax - 0.5, pen));
    check(std::fabs(sl.m.bodies[0].v(0) + v0) < 0.01*v0, fmt("elastic rebound: v = %.4f m/s", sl.m.bodies[0].v(0)));
}

static void test_latching()
{
    std::cout<<"latching and declutching"<<std::endl;
    const double mass = 10.0, k = 1000.0, T = 2.0*M_PI*std::sqrt(mass/k);
    {
        pto_settings s;
        s.control = pto_control::latching;
        s.t_hold = 0.3*T;
        s.c_latch = 1.0e6;
        s.v_eps = 1.0e-3;
        slider sl(mass, k, s);
        sl.m.bodies[0].v(0) = 0.5;
        timeint ti;
        ti.initialize(sl.m, 0.0);
        const double h = T/400.0;
        double vmax_held = 0.0, xturn = 0.0;
        int nheld = 0, holds = 0;
        bool was = false;
        for(int i=0; i<2000; ++i)
        {
            ti.step(sl.m, i*h, h);
            if(sl.pto->held())
            {
                if(!was) {++holds; xturn = sl.m.bodies[0].x(0);}
                ++nheld;
                if(nheld > 2) vmax_held = std::max(vmax_held, std::fabs(sl.m.bodies[0].v(0)));
            }
            was = sl.pto->held();
        }
        (void)xturn;
        check(holds >= 3, fmt("held at %.0f turning points in %.1f natural periods", holds, 2000.0/400.0));
        check(vmax_held < 0.02*0.5, fmt("velocity while held below %.1e m/s", vmax_held));
        check(std::fabs(nheld*h - holds*0.3*T) < holds*2.0*h, fmt("hold time %.3f s per hold (t_hold %.3f s)", nheld*h/holds, 0.3*T));
    }
    {
        pto_settings s;
        s.c = 50.0;
        s.control = pto_control::declutching;
        s.t_hold = 0.2*T;
        slider sl(mass, k, s);
        sl.m.bodies[0].v(0) = 0.5;
        timeint ti;
        ti.initialize(sl.m, 0.0);
        const double h = T/400.0;
        double Ffree = 0.0;
        for(int i=0; i<1200; ++i)
        {
            ti.step(sl.m, i*h, h);
            if(sl.pto->held()) Ffree = std::max(Ffree, std::fabs(sl.pto->force(sl.m)));
        }
        check(Ffree==0.0, "declutched: no PTO force");
    }
}

static void test_power_limit()
{
    std::cout<<"power and force limits"<<std::endl;
    pto_settings s;
    s.c = 500.0;
    s.P_max = 50.0;
    s.F_max = 300.0;
    slider sl(10.0, 1000.0, s);
    sl.m.bodies[0].x(0) = 0.5;
    timeint ti;
    ti.initialize(sl.m, 0.0);
    double Pmax = 0.0, Fmax = 0.0;
    for(int i=0; i<2000; ++i)
    {
        ti.step(sl.m, i*1.0e-3, 1.0e-3);
        Pmax = std::max(Pmax, sl.pto->power(sl.m));
        Fmax = std::max(Fmax, std::fabs(sl.pto->force(sl.m)));
    }
    check(Pmax <= 50.0 + 1.0e-9 && Pmax > 45.0, fmt("absorbed power at most P_max (max %.3f W)", Pmax));
    check(Fmax <= 300.0 + 1.0e-9, fmt("force at most F_max (max %.2f N)", Fmax));
}

static void test_rotational()
{
    std::cout<<"rotational PTO in a revolute joint"<<std::endl;
    // flap hinged at the ground (axis y), torsion spring and damper in the PTO
    model m;
    rigid_body b;
    b.mass = 1.0;
    b.inertia = Vec3(0.5, 2.0, 0.5).asDiagonal();
    m.add_body(b);
    joint &j = static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, -1, 0, Vec3::Zero(), Vec3::UnitY())));
    pto_settings s;
    s.k = 50.0;
    s.c = 1.0;
    pto_load &p = static_cast<pto_load&>(m.add_load(std::make_unique<pto_load>(j, true, s)));
    m.bodies[0].w = Vec3(0.0, 1.0, 0.0);
    timeint ti;
    ti.initialize(m, 0.0);
    const double J = 2.0, zeta = 1.0/(2.0*std::sqrt(50.0*J)), wn = std::sqrt(50.0/J), wd = wn*std::sqrt(1.0 - zeta*zeta);
    const double tend = 1.0;
    for(int i=0; i<1000; ++i)
        ti.step(m, i*1.0e-3, 1.0e-3);
    const double th = 1.0/wd*std::exp(-zeta*wn*tend)*std::sin(wd*tend);
    check(std::fabs(p.coordinate(m) - th) < 1.0e-4, fmt("rotation %.6f rad vs damped oscillator %.6f rad", p.coordinate(m), th));
}

// two-body point absorber: float (body 0) and spar (body 1), heave only
static const double rho = 1025.0, grav = 9.81;
static double Bsyn(double w, double b, double c) {const double x = w/c; return b*x*x*std::exp(-x*x);}
static double Ksyn(double t, double b, double c) {return b*c/(2.0*std::sqrt(M_PI))*(1.0 - 0.5*c*c*t*t)*std::exp(-0.25*c*c*t*t);}
static double Asyn(double w, double Ai, double b, double c)
{
    const int n = 20000;
    const double T = 12.0/c, h = T/n;
    double s = 0.0;
    for(int i=0; i<=n; ++i)
    {
        const double t = i*h;
        s += (i==0 || i==n ? 1.0 : (i%2 ? 4.0 : 2.0))*Ksyn(t, b, c)*std::sin(w*t);
    }
    return Ai - s*h/3.0/w;
}

static void test_point_absorber()
{
    std::cout<<"two-body point absorber in regular waves against the frequency domain"<<std::endl;
    const double m1 = 7.0e5, m2 = 8.5e5, C1 = 3.0e6, C2 = 3.0e5;
    const double Ai1 = 5.0e5, Ai2 = 9.0e5, b1 = 4.0e5, b2 = 1.0e5, c1 = 1.2, c2 = 1.0;
    const int n = 160;
    auto db = std::make_shared<hydro_db>();
    db->resize(2, n, 1);
    db->rho = rho;
    db->g = grav;
    db->body[0].name = "float";
    db->body[1].name = "spar";
    db->body[1].ref = Vec3(0.0, 0.0, -20.0);
    db->C(2,2) = C1;
    db->C(8,8) = C2;
    db->Ainf(2,2) = Ai1;
    db->Ainf(8,8) = Ai2;
    db->has_ainf = true;
    for(int k=0; k<n; ++k)
    {
        const double w = 0.025*(k+1);
        db->w[k] = w;
        db->A[k](2,2) = Asyn(w, Ai1, b1, c1);
        db->A[k](8,8) = Asyn(w, Ai2, b2, c2);
        db->B[k](2,2) = Bsyn(w, b1, c1);
        db->B[k](8,8) = Bsyn(w, b2, c2);
        db->X[k][0](2) = C1*std::exp(-std::pow(w/1.5, 2))*std::exp(cplx(0.0, -0.4*w));
        db->X[k][0](8) = 0.3*C2*std::exp(-std::pow(w/1.2, 2))*std::exp(cplx(0.0, -0.6*w));
    }
    const double c = 1.2e6;     // PTO damping
    for(double w : {0.6, 1.0})
    {
        const double amp = 0.5, T = 2.0*M_PI/w;
        auto waves = std::make_shared<wave_field>(-1.0, grav);
        waves->add_component(amp, w);
        waves->set_ramp(3.0*T);
        model m;
        rigid_body bf, bs;
        bf.name = "float"; bf.mass = m1; bf.inertia = Vec3(2e7, 2e7, 3e7).asDiagonal();
        bs.name = "spar"; bs.mass = m2; bs.inertia = Vec3(9e7, 9e7, 1e7).asDiagonal(); bs.x = Vec3(0.0, 0.0, -20.0);
        m.add_body(bf);
        m.add_body(bs);
        m.add_load(std::make_unique<gravity_load>());
        hydro_settings hs;
        hs.t_memory = 12.0;
        hs.dt_memory = 0.01;
        m.add_load(std::make_unique<hydro_load>(m, db, std::vector<int>{0, 1}, waves, hs));
        m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 1, Vec3(0,0,-20), Vec3::UnitZ()));
        joint &jp = static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, 1, 0, Vec3::Zero(), Vec3::UnitZ())));
        pto_settings ps;
        ps.c = c;
        pto_load &pto = static_cast<pto_load&>(m.add_load(std::make_unique<pto_load>(jp, false, ps)));
        timeint ti;
        ti.initialize(m, 0.0);
        const double h = T/150.0;
        std::vector<double> tt, z1, z2;
        double Ea = 0.0, ta = 0.0;
        const int nsteps = (int)std::lround(40.0*T/h);
        for(int i=0; i<nsteps; ++i)
        {
            ti.step(m, i*h, h);
            const double t = (i+1)*h;
            tt.push_back(t);
            z1.push_back(m.bodies[0].x(2));
            z2.push_back(m.bodies[1].x(2) + 20.0);
            if(t > 30.0*T) {Ea += pto.power(m)*h; ta += h;}
        }
        auto fit = [&](const std::vector<double> &z)
        {
            double cc = 0, ss = 0, cs = 0, zc = 0, zs = 0;
            for(size_t i=0; i<tt.size(); ++i)
            {
                if(tt[i] < 30.0*T) continue;
                const double co = std::cos(w*tt[i]), si = std::sin(w*tt[i]);
                cc += co*co; ss += si*si; cs += co*si; zc += z[i]*co; zs += z[i]*si;
            }
            const double det = cc*ss - cs*cs;
            return cplx((zc*ss - zs*cs)/det, -(zs*cc - zc*cs)/det);
        };
        const cplx Z1 = fit(z1), Z2 = fit(z2);
        const MatX A = db->added_mass(w), B = db->damping(w);
        const VecXc X = db->excitation(w, 0.0);
        Eigen::Matrix2cd K;
        K(0,0) = C1 - w*w*(m1 + A(2,2)) + cplx(0.0, w*(B(2,2) + c));
        K(1,1) = C2 - w*w*(m2 + A(8,8)) + cplx(0.0, w*(B(8,8) + c));
        K(0,1) = K(1,0) = cplx(0.0, -w*c);
        const Eigen::Vector2cd Z = K.inverse()*Eigen::Vector2cd(amp*X(2), amp*X(8));
        const cplx R = Z1 - Z2, Rf = Z(0) - Z(1);
        const double er = std::abs(R - Rf)/std::abs(Rf);
        const double P = Ea/ta, Pf = 0.5*c*w*w*std::norm(Rf);
        check(er < 5.0e-3, fmt("w = %.1f rad/s: relative motion |Z1 - Z2| %.4f m vs %.4f m", w, std::abs(R), std::abs(Rf)) + fmt(" (rel. error %.1e)", er));
        check(std::fabs(P - Pf) < 1.0e-2*Pf, fmt("w = %.1f rad/s: mean absorbed power %.2f kW vs 1/2 c w^2 |Z1 - Z2|^2 = %.2f kW", w, P/1e3, Pf/1e3));
    }
}


// --- PTO element library ---

static void test_elements()
{
    std::cout<<"PTO elements: Jacobians, composite, limits"<<std::endl;
    // each element: dF/dq and dF/dqd against central differences at states inside and outside the stops
    std::vector<std::unique_ptr<pto_base>> el;
    el.push_back(std::make_unique<pto_damper>(300.0));
    el.push_back(std::make_unique<pto_springdamper>(-200.0, 150.0, 0.1, 40.0));
    el.push_back(std::make_unique<pto_coulomb>(80.0, 0.05));
    el.push_back(std::make_unique<pto_endstop>(-0.3, 0.3, 1.0e5, 2.0e3));
    auto latch = std::make_unique<pto_latch>(5.0e4, 1.0e3);
    latch->set_param("q0", 0.05);
    latch->set_param("latch", 1.0);
    latch->set_param("t_release", 2.0);
    el.push_back(std::move(latch));
    double err = 0.0;
    for(auto &e : el)
        for(double q : {-0.4, -0.1, 0.12, 0.35})
            for(double qd : {-0.7, 0.02, 0.5})
            {
                pto_state s{q, qd, 1.0};
                const pto_output o = e->force(s);
                const double h = 1.0e-6;
                pto_state a = s, b = s;
                a.q += h; b.q -= h;
                const double dq = (e->force(a).F - e->force(b).F)/(2.0*h);
                a = s; b = s;
                a.qd += h; b.qd -= h;
                const double dqd = (e->force(a).F - e->force(b).F)/(2.0*h);
                const double sc = 1.0 + std::fabs(o.dFdq) + std::fabs(o.dFdqd);
                err = std::max(err, std::max(std::fabs(dq - o.dFdq), std::fabs(dqd - o.dFdqd))/sc);
            }
    check(err < 1.0e-6, fmt("element Jacobians against finite differences (max rel. error %.1e)", err));
    {
        pto_endstop st(-0.3, 0.3, 1.0e5, 2.0e3);
        const double Fin = st.force(pto_state{0.31, 0.5, 0.0}).F, Fout = st.force(pto_state{0.301, -1.0, 0.0}).F;
        check(std::fabs(Fin + 1.0e5*0.01 + 2.0e3*0.5) < 1.0e-6 && Fout == 0.0, fmt("end stop pushes only (%.0f N into the stop, %.0f N leaving it)", Fin, Fout));
    }
    {
        pto_latch l(5.0e4, 1.0e3);
        l.set_param("q0", 0.05); l.set_param("latch", 1.0); l.set_param("t_release", 2.0);
        const double F1 = l.force(pto_state{0.06, 0.0, 1.99}).F, F2 = l.force(pto_state{0.06, 0.0, 2.0}).F;
        check(std::fabs(F1 + 500.0) < 1.0e-9 && F2 == 0.0, fmt("latch holds q0 until t_release (%.1f N, then %.1f N)", F1, F2));
    }
    // composite: generator (damper + spring-damper) limited, end stops unchanged
    pto_composite c;
    c.add(std::make_unique<pto_damper>(1000.0));
    c.add(std::make_unique<pto_springdamper>(500.0, 0.0));
    c.add(std::make_unique<pto_endstop>(-0.3, 0.3, 1.0e5, 0.0));
    pto_output o = c.force(pto_state{0.2, 0.4, 0.0});
    check(std::fabs(o.F + 500.0) < 1.0e-12 && std::fabs(o.P - 200.0) < 1.0e-12 && c.sat_last == 0,
          fmt("composite: F = %.1f N, P = %.1f W (free)", o.F, o.P));
    c.Fmax = 300.0;
    o = c.force(pto_state{0.4, 0.4, 0.0});
    check(std::fabs(o.F - (-300.0 - 1.0e4)) < 1.0e-9 && std::fabs(o.P - 120.0) < 1.0e-12 && c.sat_last == 1,
          fmt("force limit on the generator only: F = %.1f N, P = %.1f W", o.F, o.P));
    c.Fmax = 0.0;
    c.Pmax = 50.0;
    o = c.force(pto_state{0.0, 0.5, 0.0});
    check(std::fabs(o.P - 50.0) < 1.0e-12 && c.sat_last == 2 && std::fabs(o.dFdqd - 50.0/0.25) < 1.0e-9,
          fmt("power limit: P = %.1f W, dF/dqd = -F/qd = %.1f N s/m", o.P, o.dFdqd));
    c.Pmax = 0.0;
    c.t_off0 = 1.0; c.t_off1 = 2.0;
    const double Fa = c.force(pto_state{0.0, 0.5, 0.999}).F, Fb = c.force(pto_state{0.0, 0.5, 1.0}).F, Fc = c.force(pto_state{0.0, 0.5, 2.0}).F;
    check(Fa == -500.0 && Fb == 0.0 && Fc == -500.0, "declutching window [t_off0, t_off1): generator off");
    check(c.set_param("damper", "B", 10.0) && !c.set_param("coulomb", "Fc", 1.0) && std::fabs(c.force(pto_state{0.0, 0.5, 0.0}).F + 5.0) < 1.0e-12,
          "set_param on the elements of one type");
}

static void test_controllers()
{
    std::cout<<"PTO controllers: zero crossings, latching, declutching, B_opt"<<std::endl;
    // zero crossing of qd = cos(t) at pi/2, sampled with dt = 0.3: linear interpolation
    pto_zero_crossing zc;
    double tc = -1.0;
    for(int i=0; i<10; ++i)
    {
        const double t = 0.3*i;
        if(zc.detect(t, std::cos(t))) {tc = zc.t_cross; break;}
    }
    const double t1 = 1.5, t2 = 1.8, tl = t1 + 0.3*std::cos(t1)/(std::cos(t1) - std::cos(t2));
    check(std::fabs(tc - tl) < 1.0e-12 && std::fabs(tc - M_PI/2.0) < 0.01, fmt("zero crossing at %.4f s (pi/2 = %.4f s)", tc, M_PI/2.0));
    // latching: the latch element is set at the turn and released at t_cross + t_hold
    {
        pto_composite c;
        c.add(std::make_unique<pto_damper>(100.0));
        auto *l = new pto_latch(1.0e6, 1.0e3);
        c.add(std::unique_ptr<pto_base>(l));
        pto_ctrl_latching ctrl(0.5);
        double t_on = -1.0, t_off = -1.0;
        for(int i=0; i<80; ++i)
        {
            const double t = 0.05*i;
            ctrl.update(pto_ctrl_in{t, std::sin(t), std::cos(t)}, c);
            if(ctrl.state() > 0.5 && t_on < 0.0) t_on = t;
            if(t_on > 0.0 && ctrl.state() < 0.5 && t_off < 0.0) t_off = t;
        }
        check(t_on > M_PI/2.0 && t_on < M_PI/2.0 + 0.05 && std::fabs(l->q0 - std::sin(t_on)) < 1.0e-12 && std::fabs(l->t_release - (M_PI/2.0 + 0.5)) < 2.0e-3,
              fmt("latched at %.3f s, q0 = %.4f, released at the stage time %.4f s", t_on, l->q0, l->t_release));
        check(t_off >= l->t_release && t_off < l->t_release + 0.05, fmt("latching controller off at %.3f s", t_off));
    }
    // declutching: generator window from the step of the turn to t_cross + t_off
    {
        pto_composite c;
        c.add(std::make_unique<pto_damper>(100.0));
        pto_ctrl_declutching ctrl(0.4);
        double t_on = -1.0;
        for(int i=0; i<40; ++i)
        {
            const double t = 0.05*i;
            ctrl.update(pto_ctrl_in{t, std::sin(t), std::cos(t)}, c);
            if(ctrl.state() > 0.5 && t_on < 0.0) t_on = t;
        }
        check(t_on > M_PI/2.0 && t_on < M_PI/2.0 + 0.05 && c.t_off0 == t_on && std::fabs(c.t_off1 - (M_PI/2.0 + 0.4)) < 2.0e-3,
              fmt("declutched at %.3f s until %.4f s", t_on, c.t_off1));
    }
    // hydrodynamic table (unsorted, comma separated, header) and the optimal passive damping
    const std::string file = "coral_pto_test_hydro.dat";
    {
        std::ofstream o(file);
        o<<"# w A B\n2.0, 1500, 900\n0.5 1000 100\n1.0\t1200\t400\n";
    }
    pto_hydro_table tab;
    const bool ok = tab.read(file);
    std::remove(file.c_str());
    check(ok && tab.size() == 3 && std::fabs(tab.A(1.5) - 1350.0) < 1.0e-9 && std::fabs(tab.B(0.75) - 250.0) < 1.0e-9 && tab.B(3.0) == 900.0,
          fmt("hydrodynamic table: A(1.5) = %.1f, B(0.75) = %.1f, end value B(3) = %.0f", tab.A(1.5), tab.B(0.75), tab.B(3.0)));
    const double w = 1.0, m = 5000.0, C = 2.0e4;
    const double Bo = pto_bopt(w, m, tab.A(w), tab.B(w), C);
    // B_opt maximises the mean power 1/2 B w^2 |X|^2/|C - w^2 (m+A) + i w (B_rad + B)|^2 at constant excitation
    auto P = [&](double B) {const cplx Z = C - w*w*(m + tab.A(w)) + cplx(0.0, w*(tab.B(w) + B)); return B/std::norm(Z);};
    check(P(Bo) > P(0.99*Bo) && P(Bo) > P(1.01*Bo), fmt("B_opt = %.1f N s/m maximises the absorbed power", Bo));
}

static void test_joint()
{
    std::cout<<"PTO joint of one body to the ground"<<std::endl;
    const Vec3 axis = Vec3(0.3, -0.2, 1.0), c0 = Vec3(1.0, 2.0, -0.5), xatt = Vec3(1.4, 1.7, 0.2);
    const Mat3 R0 = exp_so3(Vec3(0.1, -0.2, 0.3));
    pto_joint_prismatic j;
    j.initialize(axis, xatt, c0, R0);
    const pto_state s0 = j.state(c0, R0, Vec3::Zero(), Vec3::Zero(), 0.0);
    // moved pose: q against the attachment point, qd against finite differences of q along the motion
    const Vec3 v(0.2, 0.5, -0.3), w(0.4, -0.1, 0.25);
    auto pose = [&](double t, Vec3 &c, Mat3 &R) {c = c0 + t*v; R = exp_so3(t*w)*R0;};
    const double t = 0.7, h = 1.0e-6;
    Vec3 c, ca, cb; Mat3 R, Ra, Rb;
    pose(t, c, R); pose(t + h, ca, Ra); pose(t - h, cb, Rb);
    const pto_state s = j.state(c, R, v, w, t);
    const double qa = j.state(ca, Ra, v, w, t).q, qb = j.state(cb, Rb, v, w, t).q;
    const Vec3 a = axis.normalized();
    const double qref = a.dot(c + R*R0.transpose()*(xatt - c0) - xatt);
    j.state(c, R, v, w, t);
    Eigen::Matrix<double,6,1> nu;
    nu << v, w;
    check(std::fabs(s0.q) < 1.0e-14 && std::fabs(s.q - qref) < 1.0e-12, fmt("q = %.6f m (0 at the start)", s.q));
    check(std::fabs(s.qd - (qa - qb)/(2.0*h)) < 1.0e-7 && std::fabs(s.qd - j.J.dot(nu)) < 1.0e-14, fmt("qd = J^T [v; w] = %.6f m/s = dq/dt", s.qd));
}

static void test_latch_spring()
{
    std::cout<<"latching with a stiff latch spring in the joint PTO"<<std::endl;
    const double mass = 10.0, k = 1000.0, T = 2.0*M_PI*std::sqrt(mass/k);
    pto_settings s;
    s.c = 5.0;
    s.control = pto_control::latching;
    s.t_hold = 0.25*T;
    s.c_latch = 2.0e3;
    s.k_latch = 1.0e6;
    slider sl(mass, k, s);
    sl.m.bodies[0].v(0) = 0.5;
    timeint ti;
    ti.initialize(sl.m, 0.0);
    const double h = T/400.0;
    double dq = 0.0, q_hold = 0.0;
    int holds = 0, nheld = 0;
    bool was = false;
    for(int i=0; i<1600; ++i)
    {
        ti.step(sl.m, i*h, h);
        const bool held = sl.pto->held();
        if(held && !was) {++holds; q_hold = sl.m.bodies[0].x(0); nheld = 0;}
        if(held && ++nheld > 5) dq = std::max(dq, std::fabs(sl.m.bodies[0].x(0) - q_hold));
        was = held;
    }
    // the latch spring holds the position of the turn against the spring force k q
    check(holds >= 3 && dq < 2.0*k*0.08/1.0e6, fmt("held at %.0f turns, drift from the latched position %.1e m", holds, dq));
}

int main()
{
    test_damper();
    test_coulomb();
    test_end_stop();
    test_latching();
    test_power_limit();
    test_rotational();
    test_point_absorber();
    test_elements();
    test_controllers();
    test_joint();
    test_latch_spring();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
