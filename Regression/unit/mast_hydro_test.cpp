// Standalone tests of the REEF3D::MAST engineering hydrodynamics (src/mast/hydro): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         mast_hydro_test.cpp ../../src/mast/core/*.cpp ../../src/mast/structure/*.cpp
//         ../../src/mast/hydro/*.cpp -o mast_hydro_test
// Run:    ./mast_hydro_test
//
// A synthetic database with an analytic radiation kernel (B(w) = b (w/c)^2 exp(-(w/c)^2),
// K(t) = b c/(2 sqrt(pi)) (1 - c^2 t^2/2) exp(-c^2 t^2/4)) and A(w) from the Kramers-Kronig relation.
// Tests: waves (dispersion, spectrum Hs, kinematics consistent with the elevation), retardation
// function against the analytic kernel, infinite-frequency added mass (Ogilvie), state-space
// realization of K, database file round trip, WAMIT reader (nondimensional files), heave response in
// regular waves against the frequency-domain RAO (convolution and state-space), decay with
// convolution against state-space, hydrostatic pitch stiffness with the reference point away from
// the centre of gravity (buoyancy at the reference point), two bodies with hydrodynamic coupling
// against the 2-DOF frequency-domain solution; Morison: inertia force on a submerged vertical
// cylinder in deep-water waves, drag in a current, added mass of a free member (Newton with
// acceleration-dependent loads).
// Architect: Hans Bihs
#include"mast/core/mast_model.h"
#include"mast/core/mast_timeint.h"
#include"mast/core/mast_gravity.h"
#include"mast/core/mast_rotation.h"
#include"mast/hydro/mast_waves.h"
#include"mast/hydro/mast_hydro_db.h"
#include"mast/hydro/mast_radiation.h"
#include"mast/hydro/mast_hydro_load.h"
#include"mast/hydro/mast_morison.h"
#include<iostream>
#include<fstream>
#include<string>
#include<cstdio>
#include<cmath>
#include<complex>
#include<functional>
#include<memory>

using namespace mast;
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

// synthetic heave body: cylinder-like numbers
static const double rho = 1025.0, grav = 9.81;
static const double C33 = 31600.0, mass = 3220.0, Ainf33 = 1500.0, bK = 2000.0, cK = 2.0;

static double Bsyn(double w) {const double x = w/cK; return bK*x*x*std::exp(-x*x);}
static double Ksyn(double t) {return bK*cK/(2.0*std::sqrt(M_PI))*(1.0 - 0.5*cK*cK*t*t)*std::exp(-0.25*cK*cK*t*t);}
static double Asyn(double w)
{
    // A(w) = A_inf - 1/w int_0^inf K(t) sin(w t) dt (Simpson, K negligible beyond 12/c)
    const int n = 20000;
    const double T = 12.0/cK, h = T/n;
    double s = 0.0;
    for(int i=0; i<=n; ++i)
    {
        const double t = i*h;
        s += (i==0 || i==n ? 1.0 : (i%2 ? 4.0 : 2.0))*Ksyn(t)*std::sin(w*t);
    }
    return Ainf33 - s*h/3.0/w;
}
static cplx Xsyn(double w) {return C33*std::exp(-std::pow(w/2.5, 2))*std::exp(cplx(0.0, -0.3*w));}

// one body, coefficients only in heave (the other modes are constrained in the tests)
static std::shared_ptr<hydro_db> make_db(bool with_ainf)
{
    auto db = std::make_shared<hydro_db>();
    const int n = 160;
    db->resize(1, n, 1);
    db->rho = rho;
    db->g = grav;
    db->body[0].name = "buoy";
    db->body[0].ref = Vec3::Zero();
    db->body[0].volume = mass/rho;
    db->C(2,2) = C33;
    for(int k=0; k<n; ++k)
    {
        const double w = 0.05*(k+1);
        db->w[k] = w;
        db->A[k](2,2) = Asyn(w);
        db->B[k](2,2) = Bsyn(w);
        db->X[k][0](2) = Xsyn(w);
    }
    if(with_ainf)
    {
        db->Ainf(2,2) = Ainf33;
        db->has_ainf = true;
    }
    return db;
}

static void test_waves()
{
    std::cout<<"waves"<<std::endl;
    const double w = 1.3;
    check(std::fabs(wave_field::wavenumber(w, -1.0, grav) - w*w/grav) < 1.0e-14, "deep water k = w^2/g");
    const double kd = wave_field::wavenumber(0.2, 5.0, grav);
    check(std::fabs(0.2/kd - std::sqrt(grav*5.0)) < 0.01*std::sqrt(grav*5.0), fmt("shallow water c = sqrt(g d) (%.4f vs %.4f m/s)", 0.2/kd, std::sqrt(grav*5.0)));
    const double k20 = wave_field::wavenumber(w, 20.0, grav);
    check(std::fabs(grav*k20*std::tanh(k20*20.0) - w*w) < 1.0e-12, "finite depth: dispersion residual");

    wave_field js(-1.0, grav);
    js.add_spectrum(wave_spectrum::jonswap, 2.5, 9.0, 3.3, 0.0, 0.2, 2.5, 300, 7u);
    check(std::fabs(js.hs() - 2.5) < 0.02*2.5, fmt("JONSWAP: Hs of the components %.3f m (2.5 m)", js.hs()));
    wave_field js2(-1.0, grav);
    js2.add_spectrum(wave_spectrum::jonswap, 2.5, 9.0, 3.3, 0.0, 0.2, 2.5, 300, 7u);
    check(js2.components()[17].phase==js.components()[17].phase, "same seed, same phases");

    // regular wave in finite depth: w = d eta/dt at the surface, Bernoulli p/rho = g eta at z = 0
    wave_field rw(10.0, grav);
    rw.add_regular(1.0, 6.0, 0.4, 0.3);
    const Vec3 p(3.0, -2.0, 0.0);
    const double t = 2.7, dt = 1.0e-5;
    Vec3 u, du; double pr;
    rw.kinematics(p, t, u, du, pr);
    const double deta = (rw.elevation(p(0), p(1), t + dt) - rw.elevation(p(0), p(1), t - dt))/(2.0*dt);
    check(std::fabs(u(2) - deta) < 1.0e-7, fmt("vertical velocity at z = 0 = d eta/dt (%.6f vs %.6f)", u(2), deta));
    check(std::fabs(pr - grav*rw.elevation(p(0), p(1), t)) < 1.0e-12, "dynamic pressure at z = 0 = g eta");
    Vec3 u1, du1, u2, du2; double q;
    rw.kinematics(p, t + dt, u1, du1, q);
    rw.kinematics(p, t - dt, u2, du2, q);
    check(((u1 - u2)/(2.0*dt) - du).norm() < 1.0e-6, "acceleration = du/dt");
    // horizontal velocity along the direction of propagation
    check(std::fabs(u(0)*std::sin(0.4) - u(1)*std::cos(0.4)) < 1.0e-14, "horizontal velocity along beta");
}

static void test_retardation()
{
    std::cout<<"radiation kernel"<<std::endl;
    auto db = make_db(false);
    retardation K;
    K.compute(*db, 10.0, 0.01);
    double err = 0.0;
    for(int s=0; s<K.nsample(); ++s)
        err = std::max(err, std::fabs(K.sample(s)(2,2) - Ksyn(s*0.01)));
    const double kmax = Ksyn(0.0);
    check(err < 2.0e-3*kmax, fmt("K(t) from B(w) against the analytic kernel: max error %.1e of K(0)", err/kmax));
    const MatX Ai = ogilvie_ainf(*db, K);
    check(std::fabs(Ai(2,2) - Ainf33) < 0.005*Ainf33, fmt("A_inf (Ogilvie) %.1f kg (exact %.1f kg)", Ai(2,2), Ainf33));
    radiation_ss ss;
    ss.identify(K, 1.0e-5, 12);
    check(ss.pairs.size()==1 && ss.pairs[0].i==2 && ss.pairs[0].j==2, "state-space: one pair (heave-heave)");
    if(!ss.pairs.empty())
    {
        check(ss.pairs[0].error < 5.0e-3, fmt("state-space: order %.0f, fit error %.1e of max |K|", (double)ss.pairs[0].lam.size(), ss.pairs[0].error));
        bool stable = true;
        for(int r=0; r<ss.pairs[0].lam.size(); ++r)
            stable = stable && ss.pairs[0].lam(r).real() < 0.0;
        check(stable, "state-space: all poles stable");
    }
}

static void test_db_files()
{
    std::cout<<"database files"<<std::endl;
    auto db = make_db(true);
    db->write("/tmp/mast_hydro_test.mhd");
    hydro_db rd;
    rd.read("/tmp/mast_hydro_test.mhd");
    double e = std::fabs(rd.C(2,2) - db->C(2,2)) + std::fabs(rd.Ainf(2,2) - db->Ainf(2,2)) + std::fabs(rd.body[0].volume - db->body[0].volume);
    for(int k=0; k<db->nfreq(); ++k)
        e += std::fabs(rd.w[k] - db->w[k]) + (rd.A[k] - db->A[k]).norm() + (rd.B[k] - db->B[k]).norm() + (rd.X[k][0] - db->X[k][0]).norm();
    check(e < 1.0e-8 && rd.body[0].name=="buoy" && rd.has_ainf, fmt("MAST format: write and read back (difference %.1e)", e));
    const double wq = 1.234;
    check(std::fabs(rd.added_mass(wq)(2,2) - Asyn(wq)) < 0.5 && std::fabs(rd.damping(wq)(2,2) - Bsyn(wq)) < 1.0
          && std::abs(rd.excitation(wq, 0.0)(2) - Xsyn(wq)) < 5.0, "interpolation between frequencies");

    // WAMIT: nondimensional files with ULEN = 2
    const double L = 2.0;
    {
        std::ofstream f1("/tmp/mast_wamit.1"), f3("/tmp/mast_wamit.3"), fh("/tmp/mast_wamit.hst");
        f1<<"  -1.0  3  3  0.5\n   0.0  3  3  0.9\n";
        f3<<std::scientific;
        for(double P : {8.0, 4.0})
        {
            const double om = 2.0*M_PI/P;
            f1<<P<<"  3  3  "<<0.6<<"  "<<0.2<<"\n"<<P<<"  5  5  "<<0.03<<"  "<<0.01<<"\n"<<P<<"  3  5  "<<0.02<<"  "<<0.004<<"\n";
            f3<<P<<"  0.0  3  1.0  -30.0  "<<std::cos(-M_PI/6.0)<<"  "<<std::sin(-M_PI/6.0)<<"\n";
            f3<<P<<"  0.0  5  0.1  0.0  0.1  0.0\n";
            (void)om;
        }
        fh<<"  3  3  3.14\n  5  5  0.8\n";
    }
    hydro_db wd;
    wd.read_wamit("/tmp/mast_wamit", 1, L, rho, grav);
    const double w8 = 2.0*M_PI/8.0;
    bool ok = wd.nfreq()==2 && std::fabs(wd.w[0] - w8) < 1.0e-12;
    ok = ok && std::fabs(wd.A[0](2,2) - 0.6*rho*std::pow(L,3)) < 1.0e-9;
    ok = ok && std::fabs(wd.B[0](2,2) - 0.2*rho*w8*std::pow(L,3)) < 1.0e-9;
    ok = ok && std::fabs(wd.A[0](4,4) - 0.03*rho*std::pow(L,5)) < 1.0e-9;
    ok = ok && std::fabs(wd.A[0](2,4) - 0.02*rho*std::pow(L,4)) < 1.0e-9;
    ok = ok && wd.has_ainf && std::fabs(wd.Ainf(2,2) - 0.5*rho*std::pow(L,3)) < 1.0e-9;
    ok = ok && std::abs(wd.X[0][0](2) - std::polar(rho*grav*L*L, -M_PI/6.0)) < 1.0e-6*rho*grav*L*L;
    ok = ok && std::abs(wd.X[0][0](4) - 0.1*rho*grav*std::pow(L,3)) < 1.0e-6;
    ok = ok && std::fabs(wd.C(2,2) - 3.14*rho*grav*L*L) < 1.0e-6 && std::fabs(wd.C(4,4) - 0.8*rho*grav*std::pow(L,4)) < 1.0e-6;
    check(ok, "WAMIT .1/.3/.hst: periods, A, B, A_inf, X and C made dimensional (ULEN, rho, g)");
}

// heave-only body in a hydrodynamic database; returns the model, waves set by the caller
struct heave_case
{
    model m;
    std::shared_ptr<wave_field> waves;
    hydro_load *hl = nullptr;
};

static void build_heave(heave_case &hc, std::shared_ptr<hydro_db> db, std::shared_ptr<wave_field> waves, bool ss, double Bpto)
{
    rigid_body b;
    b.name = "buoy";
    b.mass = mass;
    b.inertia = Vec3(5000.0, 5000.0, 3000.0).asDiagonal();
    hc.m.add_body(b);
    hc.m.add_load(std::make_unique<gravity_load>());
    hydro_settings s;
    s.state_space = ss;
    s.t_memory = 8.0;
    s.dt_memory = 0.005;
    s.linear_damping = MatX::Zero(6, 6);
    s.linear_damping(2,2) = Bpto;
    hc.hl = &static_cast<hydro_load&>(hc.m.add_load(std::make_unique<hydro_load>(hc.m, db, std::vector<int>{0}, waves, s)));
    hc.m.add_constraint(std::make_unique<joint>(hc.m, joint_type::prismatic, -1, 0, Vec3::Zero(), Vec3::UnitZ()));
    hc.waves = waves;
}

// amplitude and phase of z(t) = Re{Z exp(i w t)} over [t0, t1] by least squares
static cplx fit(const std::vector<double> &t, const std::vector<double> &z, double w, double t0)
{
    double cc = 0, ss = 0, cs = 0, zc = 0, zs = 0;
    for(size_t i=0; i<t.size(); ++i)
    {
        if(t[i] < t0) continue;
        const double c = std::cos(w*t[i]), s = std::sin(w*t[i]);
        cc += c*c; ss += s*s; cs += c*s; zc += z[i]*c; zs += z[i]*s;
    }
    const double det = cc*ss - cs*cs;
    const double a = (zc*ss - zs*cs)/det, b = (zs*cc - zc*cs)/det;    // z = a cos - b' sin: Z = a + i b'
    return cplx(a, -b);
}

static void test_rao(bool ss)
{
    std::cout<<"heave in regular waves against the frequency-domain RAO ("<<(ss ? "state-space" : "convolution")<<")"<<std::endl;
    auto db = make_db(true);
    const double Bpto = 800.0, amp = 0.5;
    for(double w : {0.8, 1.6, 2.6})
    {
        auto waves = std::make_shared<wave_field>(-1.0, grav);
        waves->add_component(amp, w, 0.0, 0.0);
        const double T = 2.0*M_PI/w;
        waves->set_ramp(3.0*T);
        heave_case hc;
        build_heave(hc, db, waves, ss, Bpto);
        timeint ti;
        const double h = T/120.0;
        ti.initialize(hc.m, 0.0);
        std::vector<double> tt, zz;
        const int n = (int)std::lround(30.0*T/h);
        for(int i=0; i<n; ++i)
        {
            if(!ti.step(hc.m, i*h, h)) {std::cout<<"  step failed"<<std::endl; break;}
            tt.push_back((i+1)*h);
            zz.push_back(hc.m.bodies[0].x(2));
        }
        const cplx Z = fit(tt, zz, w, 20.0*T);
        const cplx Zfd = amp*Xsyn(w)/(C33 - w*w*(mass + Asyn(w)) + cplx(0.0, w*(Bsyn(w) + Bpto)));
        const double e = std::abs(Z - Zfd)/std::abs(Zfd);
        check(e < 7.0e-3, fmt("w = %.1f rad/s: |Z| %.5f m vs %.5f m, phase error %.4f rad", w, std::abs(Z), std::abs(Zfd), std::arg(Z/Zfd)) + fmt(" (rel. error %.1e)", e));
    }
}

static void test_decay()
{
    std::cout<<"heave decay: convolution against state-space"<<std::endl;
    auto db = make_db(true);
    std::vector<double> z[2];
    for(int k=0; k<2; ++k)
    {
        heave_case hc;
        build_heave(hc, db, nullptr, k==1, 0.0);
        hc.m.bodies[0].x(2) = 0.2;
        timeint ti;
        ti.initialize(hc.m, 0.0);
        const double h = 0.01;
        for(int i=0; i<2000; ++i)
        {
            ti.step(hc.m, i*h, h);
            z[k].push_back(hc.m.bodies[0].x(2));
        }
    }
    double d = 0.0, zm = 0.0;
    for(size_t i=0; i<z[0].size(); ++i)
    {
        d = std::max(d, std::fabs(z[0][i] - z[1][i]));
        zm = std::max(zm, std::fabs(z[0][i]));
    }
    check(d < 2.0e-3*0.2, fmt("max difference %.1e m (initial offset 0.2 m)", d));
    // logarithmic decrement of the positive peaks against B(w_n) / (2 (m + A(w_n)) w_n)
    std::vector<double> peaks;
    for(size_t i=1; i+1<z[0].size(); ++i)
        if(z[0][i] > z[0][i-1] && z[0][i] >= z[0][i+1] && z[0][i] > 0.0)
            peaks.push_back(z[0][i]);
    const double wn = std::sqrt(C33/(mass + Ainf33));
    double wn2 = wn;
    for(int it=0; it<20; ++it) wn2 = std::sqrt(C33/(mass + Asyn(wn2)));
    const double zeta_ex = Bsyn(wn2)/(2.0*(mass + Asyn(wn2))*wn2);
    const double zeta = peaks.size() > 4 ? std::log(peaks[1]/peaks[4])/(3.0*2.0*M_PI) : 0.0;
    check(std::fabs(zeta - zeta_ex) < 0.2*zeta_ex, fmt("damping ratio from the decay %.4f vs B/(2 (m+A) w_n) = %.4f", zeta, zeta_ex));
    double tcross = -1.0;
    for(size_t i=1; i<z[0].size(); ++i)
        if(z[0][i-1] > 0.0 && z[0][i] <= 0.0) {tcross = 0.01*i; break;}
    const double Tn = 2.0*M_PI/wn2;
    check(std::fabs(4.0*tcross - Tn) < 0.1*Tn, fmt("first zero crossing at %.3f s ~ T_n/4 = %.3f s", tcross, Tn/4.0));
}

static void test_pitch_stiffness()
{
    std::cout<<"hydrostatic pitch stiffness with the reference point away from the CoG"<<std::endl;
    // database body with the reference point 0.5 m above the CoG; C55 about the reference point
    auto db = std::make_shared<hydro_db>();
    db->resize(1, 2, 1);
    db->rho = rho;
    db->g = grav;
    db->w = {0.5, 1.0};
    const Vec3 G(0.0, 0.0, -1.0), O(0.0, 0.0, -0.5);
    db->body[0].ref = O;
    db->body[0].volume = mass/rho;
    const double C55 = 40000.0;
    db->C(4,4) = C55;
    model m;
    rigid_body b;
    b.mass = mass;
    b.inertia = Vec3(4000.0, 4000.0, 3000.0).asDiagonal();
    b.x = G;
    m.add_body(b);
    m.add_load(std::make_unique<gravity_load>());
    hydro_settings s;
    s.radiation = false;
    s.excitation = false;
    s.linear_damping = MatX::Zero(6, 6);
    s.linear_damping(4,4) = 20000.0;
    m.add_load(std::make_unique<hydro_load>(m, db, std::vector<int>{0}, nullptr, s));
    m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, -1, 0, G, Vec3::UnitY()));
    const double M = 500.0;
    struct moment_load : public load
    {
        double M;
        explicit moment_load(double M_) : M(M_) {}
        std::string name() const override {return "moment";}
        void eval(const model &mm, double, VecX &f) const override {mm.add_moment(f, 0, Vec3(0.0, M, 0.0));}
    };
    m.add_load(std::make_unique<moment_load>(M));
    timeint ti;
    ti.initialize(m, 0.0);
    for(int i=0; i<4000; ++i)
        ti.step(m, i*0.01, 0.01);
    const double theta = log_so3(m.bodies[0].R)(1);
    // about G: C55_G = C55_O + rho g V (z_O - z_G)
    const double C55G = C55 + rho*grav*(mass/rho)*(O(2) - G(2));
    check(std::fabs(theta - M/C55G) < 1.0e-3*M/C55G, fmt("static pitch %.6e rad = M / (C55_O + B (z_O - z_G)) = %.6e rad", theta, M/C55G));
}

static void test_two_bodies()
{
    std::cout<<"two bodies with hydrodynamic coupling against the 2-DOF frequency domain"<<std::endl;
    // database with two bodies, heave only: coupling A13, B13 (body 1 heave - body 2 heave), X for both
    const int n = 160;
    auto db = std::make_shared<hydro_db>();
    db->resize(2, n, 1);
    db->rho = rho;
    db->g = grav;
    db->body[0].ref = Vec3::Zero();
    db->body[1].ref = Vec3(0.0, 0.0, -3.0);
    db->body[0].volume = mass/rho;
    db->body[1].volume = 2.0*mass/rho;
    db->C(2,2) = C33;
    db->C(8,8) = 0.5*C33;
    db->Ainf(2,2) = Ainf33; db->Ainf(8,8) = 2.0*Ainf33; db->Ainf(2,8) = db->Ainf(8,2) = 0.3*Ainf33;
    db->has_ainf = true;
    for(int k=0; k<n; ++k)
    {
        const double w = 0.05*(k+1);
        db->w[k] = w;
        const double dA = Asyn(w) - Ainf33;
        db->A[k](2,2) = Asyn(w);
        db->A[k](8,8) = 2.0*Ainf33 + 1.5*dA;
        db->A[k](2,8) = db->A[k](8,2) = 0.3*Ainf33 + 0.4*dA;
        db->B[k](2,2) = Bsyn(w);
        db->B[k](8,8) = 1.5*Bsyn(w);
        db->B[k](2,8) = db->B[k](8,2) = 0.4*Bsyn(w);
        db->X[k][0](2) = Xsyn(w);
        db->X[k][0](8) = 0.6*Xsyn(w)*std::exp(cplx(0.0, 0.5));
    }
    const double w = 1.4, amp = 0.4, T = 2.0*M_PI/w;
    auto waves = std::make_shared<wave_field>(-1.0, grav);
    waves->add_component(amp, w);
    waves->set_ramp(3.0*T);
    model m;
    rigid_body b1, b2;
    b1.mass = mass; b1.inertia = Vec3(5000,5000,3000).asDiagonal();
    b2.mass = 2.0*mass; b2.inertia = Vec3(9000,9000,5000).asDiagonal(); b2.x = Vec3(0.0, 0.0, -3.0);
    m.add_body(b1);
    m.add_body(b2);
    m.add_load(std::make_unique<gravity_load>());
    hydro_settings s;
    s.t_memory = 8.0;
    s.dt_memory = 0.005;
    s.linear_damping = MatX::Zero(12, 12);
    s.linear_damping(2,2) = 500.0;
    s.linear_damping(8,8) = 3000.0;
    m.add_load(std::make_unique<hydro_load>(m, db, std::vector<int>{0, 1}, waves, s));
    m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 0, Vec3::Zero(), Vec3::UnitZ()));
    m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 1, Vec3(0,0,-3), Vec3::UnitZ()));
    timeint ti;
    ti.initialize(m, 0.0);
    const double h = T/120.0;
    std::vector<double> tt, z1, z2;
    for(int i=0; i<(int)std::lround(30.0*T/h); ++i)
    {
        ti.step(m, i*h, h);
        tt.push_back((i+1)*h);
        z1.push_back(m.bodies[0].x(2));
        z2.push_back(m.bodies[1].x(2) + 3.0);
    }
    const cplx Z1 = fit(tt, z1, w, 20.0*T), Z2 = fit(tt, z2, w, 20.0*T);
    Eigen::Matrix2cd Kd;
    const MatX A = db->added_mass(w), B = db->damping(w);
    Kd(0,0) = C33 - w*w*(mass + A(2,2)) + cplx(0.0, w*(B(2,2) + 500.0));
    Kd(1,1) = 0.5*C33 - w*w*(2.0*mass + A(8,8)) + cplx(0.0, w*(B(8,8) + 3000.0));
    Kd(0,1) = -w*w*A(2,8) + cplx(0.0, w*B(2,8));
    Kd(1,0) = -w*w*A(8,2) + cplx(0.0, w*B(8,2));
    const VecXc Xw = db->excitation(w, 0.0);
    const Eigen::Vector2cd Zfd = Kd.inverse()*Eigen::Vector2cd(amp*Xw(2), amp*Xw(8));
    const double e1 = std::abs(Z1 - Zfd(0))/std::abs(Zfd(0)), e2 = std::abs(Z2 - Zfd(1))/std::abs(Zfd(1));
    check(e1 < 5.0e-3 && e2 < 5.0e-3, fmt("|Z1| %.5f m (%.5f), |Z2| %.5f m (%.5f)", std::abs(Z1), std::abs(Zfd(0)), std::abs(Z2), std::abs(Zfd(1))) + fmt(", rel. errors %.1e %.1e", e1, e2));
}

static void test_morison()
{
    std::cout<<"Morison members"<<std::endl;
    model m;
    // (a) fixed vertical cylinder from z = -40 to -2 in deep-water waves, inertia only
    auto waves = std::make_shared<wave_field>(-1.0, grav);
    const double a = 0.5, w = 1.0, D = 1.0;
    waves->add_component(a, w);
    morison_member mb;
    mb.p1 = Vec3(0, 0, -40.0); mb.p2 = Vec3(0, 0, -2.0); mb.D = D; mb.Cd = 0.0; mb.Ca = 1.0; mb.Cm = 2.0; mb.n = 400;
    morison_load fixed(m, -1, {mb}, waves, rho);
    double Fmax = 0.0;
    for(int i=0; i<400; ++i)
    {
        Vec3 F, M;
        fixed.total(m, i*2.0*M_PI/w/400.0, F, M);
        Fmax = std::max(Fmax, std::fabs(F(0)));
    }
    const double k = w*w/grav, A = 0.25*M_PI*D*D;
    const double Fex = rho*2.0*A*a*w*w*(std::exp(-2.0*k) - std::exp(-40.0*k))/k;
    check(std::fabs(Fmax - Fex) < 1.0e-3*Fex, fmt("inertia force amplitude %.1f N vs %.1f N (deep water)", Fmax, Fex));
    // (b) drag in a uniform current
    morison_member md;
    md.p1 = Vec3(0, 0, -10.0); md.p2 = Vec3(0, 0, -1.0); md.D = 0.5; md.Cd = 1.2; md.n = 10;
    const double U = 1.5;
    morison_load drag(m, -1, {md}, nullptr, rho, Vec3(U, 0, 0));
    Vec3 F, M;
    drag.total(m, 0.0, F, M);
    const double Fd = 0.5*rho*1.2*0.5*9.0*U*U;
    check(std::fabs(F(0) - Fd) < 1.0e-9*Fd && std::fabs(F(2)) < 1.0e-12, fmt("drag in a current %.2f N = 1/2 rho Cd D L U^2 = %.2f N", F(0), Fd));
    // (c) free member pushed sideways in still water: a = F0 / (m + rho Ca V)
    model mf;
    rigid_body b;
    b.mass = 500.0;
    b.inertia = Vec3(1000, 1000, 100).asDiagonal();
    b.x = Vec3(0, 0, -6.0);
    mf.add_body(b);
    morison_member mm;
    mm.p1 = Vec3(0, 0, -10.0); mm.p2 = Vec3(0, 0, -2.0); mm.D = 0.6; mm.Cd = 0.0; mm.Ca = 1.0; mm.n = 8;
    mf.add_load(std::make_unique<morison_load>(mf, 0, std::vector<morison_member>{mm}, nullptr, rho));
    const double F0 = 1000.0;
    struct push : public load
    {
        double F0;
        explicit push(double f) : F0(f) {}
        std::string name() const override {return "push";}
        void eval(const model &mm2, double, VecX &f) const override {mm2.add_force(f, 0, mm2.bodies[0].x, Vec3(F0, 0, 0));}
    };
    mf.add_load(std::make_unique<push>(F0));
    timeint ti;
    ti.initialize(mf, 0.0);
    for(int i=0; i<100; ++i)
        ti.step(mf, i*0.01, 0.01);
    const double ma = rho*1.0*0.25*M_PI*0.36*8.0;
    const double xex = 0.5*F0/(500.0 + ma);
    check(std::fabs(mf.bodies[0].x(0) - xex) < 1.0e-9, fmt("added mass of the member: x(1 s) = %.6f m = F0/(2 (m + rho Ca V)) = %.6f m", mf.bodies[0].x(0), xex));
}

int main()
{
    test_waves();
    test_retardation();
    test_db_files();
    test_rao(false);
    test_rao(true);
    test_decay();
    test_pitch_stiffness();
    test_two_bodies();
    test_morison();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
