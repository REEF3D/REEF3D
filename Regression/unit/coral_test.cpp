// Standalone tests of REEF3D::CORAL (src/coral/): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         coral_test.cpp ../../src/coral/core/*.cpp ../../src/coral/structure/*.cpp -o coral_test
// Run:    ./coral_test
//
// Tests: name and version, core types; rotations (exp/log, tangent operator, Euler angles);
// generalized-alpha integrator: free fall (exact), linear oscillator (period, second order, no
// numerical damping for rho_inf = 1, damping for rho_inf < 1), torque-free asymmetric top (angular
// momentum and energy, second order); joints: pendulum on a revolute joint (against the ODE, second
// order, constraint), prismatic joints in series (lateral motion shared, axial motion free, joint force),
// fixed joint (two bodies fall as one), heavy symmetric top on a spherical joint (energy, L_z),
// cylindrical and planar joints (constraint kept, free motions); elements: tether (tension-only
// spring, energy through slack phases), linear 6x6 stiffness; modified Newton (Jacobian reuse).
// Architect: Hans Bihs
#include"coral/core/coral_types.h"
#include"coral/core/coral_version.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_gravity.h"
#include"coral/structure/coral_spring.h"
#include"coral/structure/coral_linear6.h"
#include<iostream>
#include<string>
#include<cstring>
#include<cstdio>
#include<cmath>
#include<functional>
#include<memory>

using namespace coral;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}
static std::string fmt(const char *f, double a, double b=0.0, double c=0.0)
{
    char s[256];
    std::snprintf(s, sizeof(s), f, a, b, c);
    return s;
}

// constant or state-dependent external load for the tests
class test_load : public load
{
public:
    explicit test_load(std::function<void(const model&, double, VecX&)> fn) : fn_(std::move(fn)) {}
    std::string name() const override {return "test";}
    void eval(const model &m, double t, VecX &f) const override {fn_(m, t, f);}
private:
    std::function<void(const model&, double, VecX&)> fn_;
};

static rigid_body make_body(double mass, const Vec3 &J, const Vec3 &x, const std::string &name="body")
{
    rigid_body b;
    b.name = name;
    b.mass = mass;
    b.inertia = J.asDiagonal();
    b.x = x;
    return b;
}

static void run(model &m, timeint &ti, double t0, double t1, double h,
                const std::function<void(double)> &each = nullptr)
{
    ti.initialize(m, t0);
    const int n = (int)std::lround((t1 - t0)/h);
    for(int i=0; i<n; ++i)
    {
        const double t = t0 + i*h;
        if(!ti.step(m, t, h))
        {
            std::cout<<"  step failed at t = "<<t<<std::endl;
            return;
        }
        if(each) each(t + h);
    }
}

static void test_version()
{
    std::cout<<"version"<<std::endl;
    check(std::string(coral::name())=="REEF3D::CORAL", "name is REEF3D::CORAL");
    check(std::strlen(coral::version())>0, "version set");
    const std::string b = coral::banner();
    check(b.find(coral::name())==0 && b.find(coral::version())!=std::string::npos
          && b.find(coral::long_name())!=std::string::npos, "banner = name version (long name): "+b);
}

static void test_types()
{
    std::cout<<"core types"<<std::endl;
    check(Vec6::RowsAtCompileTime==6 && Mat6::ColsAtCompileTime==6, "Vec6/Mat6 are 6x1/6x6");
    Mat3 R = Quat(Eigen::AngleAxisd(0.3, Vec3::UnitZ())).toRotationMatrix();
    check((R.transpose()*R - Mat3::Identity()).norm()<1.0e-14, "Quat -> Mat3 orthonormal");
    check(sizeof(coral::index)==sizeof(std::ptrdiff_t), "index is ptrdiff_t");
}

static void test_rotation()
{
    std::cout<<"rotations"<<std::endl;
    double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0;
    const Vec3 dirs[4] = {Vec3(1,2,3).normalized(), Vec3(-0.3,0.1,0.9).normalized(), Vec3::UnitX(), Vec3(0.5,-0.5,0.7).normalized()};
    const double angles[6] = {0.0, 1.0e-9, 1.0e-4, 0.7, 2.5, 3.1415};
    for(const Vec3 &d : dirs)
        for(double a : angles)
        {
            const Vec3 psi = a*d;
            const Mat3 R = exp_so3(psi);
            e1 = std::max(e1, (R - Eigen::AngleAxisd(a, d).toRotationMatrix()).norm());
            e2 = std::max(e2, (log_so3(R) - psi).norm());
            // tangent operator: R(psi + eps dpsi) = R(psi) exp(eps T dpsi) + O(eps^2)
            const Vec3 dpsi(0.3, -0.2, 0.5);
            const double eps = 1.0e-6;
            const Vec3 fd = log_so3(R.transpose()*exp_so3(psi + eps*dpsi))/eps;
            if(a < 3.0)
                e3 = std::max(e3, (fd - tangent_so3(psi)*dpsi).norm());
            const Vec3 rpy(0.3*a, -0.2*a, 0.4*a);
            e4 = std::max(e4, (euler_from_matrix(matrix_from_euler(rpy)) - rpy).norm());
        }
    check(e1 < 1.0e-14, fmt("exp = Rodrigues (err %.1e)", e1));
    check(e2 < 1.0e-12, fmt("log(exp(psi)) = psi up to pi (err %.1e)", e2));
    check(e3 < 1.0e-5, fmt("tangent operator against finite differences (err %.1e)", e3));
    check(e4 < 1.0e-13, fmt("Euler angles round trip (err %.1e)", e4));
}

static void test_free_fall()
{
    std::cout<<"free fall (constant acceleration is exact)"<<std::endl;
    for(double rho : {1.0, 0.8, 0.3})
    {
        model m;
        rigid_body b = make_body(2.0, Vec3(1,1,1), Vec3(0,0,10));
        b.v = Vec3(1.0, 0.0, 2.0);
        m.add_body(b);
        m.add_load(std::make_unique<gravity_load>());
        timeint_settings s;
        s.rho_inf = rho;
        timeint ti(s);
        run(m, ti, 0.0, 1.0, 0.05);
        const Vec3 ex = Vec3(1.0, 0.0, 10.0 + 2.0 - 0.5*9.81);
        const double err = (m.bodies[0].x - ex).norm() + (m.bodies[0].v - Vec3(1.0, 0.0, 2.0 - 9.81)).norm();
        check(err < 1.0e-10, fmt("rho_inf %.1f: position and velocity exact (err %.1e)", rho, err));
    }
}

// mass on a linear spring (linear6), x(t) = A cos(omega t); amp = amplitude from the energy at tend
static double oscillator(double rho, double h, double tend, double &amp)
{
    model m;
    const double mass = 3.0, k = 12.0, A = 0.1;
    m.add_body(make_body(mass, Vec3(1,1,1), Vec3::Zero()));
    Mat6 K = Mat6::Zero();
    K(0,0) = k;
    m.add_load(std::make_unique<linear6>(m, 0, Vec3::Zero(), K, Mat6::Zero()));
    m.bodies[0].x = Vec3(A,0,0);
    timeint_settings s;
    s.rho_inf = rho;
    timeint ti(s);
    run(m, ti, 0.0, tend, h);
    amp = std::sqrt(2.0*(m.kinetic_energy() + m.potential_energy())/k);
    const double w = std::sqrt(k/mass);
    return std::fabs(m.bodies[0].x(0) - A*std::cos(w*tend));
}

static void test_oscillator()
{
    std::cout<<"linear oscillator"<<std::endl;
    const double T = M_PI;    // omega = 2
    double amp;
    const double e1 = oscillator(1.0, T/50.0, T, amp);
    const double e2 = oscillator(1.0, T/100.0, T, amp);
    const double e3 = oscillator(1.0, T/200.0, T, amp);
    const double p1 = std::log2(e1/e2), p2 = std::log2(e2/e3);
    // after a full period the phase errors of the trapezoidal-type scheme cancel to higher order
    check(p1 > 1.9 && p2 > 1.9, fmt("rho_inf 1: at least second order (orders %.2f %.2f)", p1, p2));
    const double f1 = oscillator(1.0, T/50.0, 0.3*T, amp), f2 = oscillator(1.0, T/100.0, 0.3*T, amp);
    check(std::fabs(std::log2(f1/f2) - 2.0) < 0.1, fmt("rho_inf 1: second order at t = 0.3 T (order %.2f)", std::log2(f1/f2)));
    check(e3/0.1 < 1.0e-3, fmt("rho_inf 1: error after one period with T/200: %.1e of the amplitude", e3/0.1));
    oscillator(1.0, T/20.0, 50.0*T, amp);
    check(std::fabs(amp - 0.1) < 1.0e-12, fmt("rho_inf 1: no numerical damping over 50 periods (amplitude %.13f)", amp));
    double amp2;
    oscillator(0.5, T/20.0, 50.0*T, amp2);
    check(amp2 < 0.0999 && amp2 > 0.07, fmt("rho_inf 0.5: numerical damping at T/20 (amplitude %.5f after 50 periods)", amp2));
    oscillator(0.5, T/200.0, 50.0*T, amp2);
    check(std::fabs(amp2 - 0.1) < 1.0e-4, fmt("rho_inf 0.5: damping vanishes with small steps (amplitude %.6f at T/200)", amp2));
}

static void test_top()
{
    std::cout<<"torque-free asymmetric top"<<std::endl;
    auto top = [](double h, double tend, Vec3 &L, double &E, Mat3 &R)
    {
        model m;
        rigid_body b = make_body(1.0, Vec3(1.0, 2.0, 3.0), Vec3::Zero());
        b.w = Vec3(3.0, 0.2, 0.1);
        m.add_body(b);
        timeint_settings s;
        s.rho_inf = 1.0;
        timeint ti(s);
        run(m, ti, 0.0, tend, h);
        L = m.bodies[0].angular_momentum();
        E = m.bodies[0].kinetic_energy();
        R = m.bodies[0].R;
    };
    rigid_body b0 = make_body(1.0, Vec3(1.0, 2.0, 3.0), Vec3::Zero());
    b0.w = Vec3(3.0, 0.2, 0.1);
    const Vec3 L0 = b0.angular_momentum();
    const double E0 = b0.kinetic_energy();
    Vec3 L1, L2, L3, Lr; double E1, E2, E3, Er; Mat3 R1, R2, R3, Rr;
    top(0.02, 2.0, L1, E1, R1);
    top(0.01, 2.0, L2, E2, R2);
    top(0.005, 2.0, L3, E3, R3);
    top(0.0005, 2.0, Lr, Er, Rr);
    check((L3 - L0).norm()/L0.norm() < 1.0e-5, fmt("angular momentum (global) kept: rel. error %.1e", (L3 - L0).norm()/L0.norm()));
    check(std::fabs(E3 - E0)/E0 < 1.0e-5, fmt("kinetic energy kept: rel. error %.1e", std::fabs(E3 - E0)/E0));
    const double p = std::log2((R1 - Rr).norm()/(R2 - Rr).norm()), q = std::log2((R2 - Rr).norm()/(R3 - Rr).norm());
    check(p > 1.8 && q > 1.8, fmt("attitude second order (orders %.2f %.2f)", p, q));
}

// pendulum: revolute joint at the origin, axis y, bob at distance L; angle from the horizontal x axis
static void test_pendulum()
{
    std::cout<<"pendulum on a revolute joint"<<std::endl;
    const double L = 1.5, mass = 2.0, r = 0.1, g = 9.81;
    const double Jyy = 0.4*mass*r*r;
    // reference: (m L^2 + J) th'' = -m g L cos th, th(0) = 0 (horizontal), RK4 with small steps
    auto ref = [&](double tend)
    {
        double th = 0.0, om = 0.0;
        const double hh = 1.0e-5;
        const int n = (int)std::lround(tend/hh);
        auto acc = [&](double t_){return -mass*g*L*std::cos(t_)/(mass*L*L + Jyy);};
        for(int i=0; i<n; ++i)
        {
            const double k1t = om, k1o = acc(th);
            const double k2t = om + 0.5*hh*k1o, k2o = acc(th + 0.5*hh*k1t);
            const double k3t = om + 0.5*hh*k2o, k3o = acc(th + 0.5*hh*k2t);
            const double k4t = om + hh*k3o, k4o = acc(th + hh*k3t);
            th += hh/6.0*(k1t + 2.0*k2t + 2.0*k3t + k4t);
            om += hh/6.0*(k1o + 2.0*k2o + 2.0*k3o + k4o);
        }
        return th;
    };
    const double tend = 1.2;
    const double thref = ref(tend);
    auto pend = [&](double h, double &cmax, double &dE)
    {
        model m;
        m.add_body(make_body(mass, Vec3(Jyy, Jyy, Jyy), Vec3(L, 0, 0)));
        m.add_load(std::make_unique<gravity_load>());
        m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, -1, 0, Vec3::Zero(), Vec3::UnitY()));
        timeint_settings s;
        s.rho_inf = 1.0;
        timeint ti(s);
        const double E0 = m.kinetic_energy() + m.potential_energy();
        run(m, ti, 0.0, tend, h);
        cmax = ti.stats.max_constraint;
        dE = std::fabs(m.kinetic_energy() + m.potential_energy() - E0);
        // th measured downwards from x: z = -L sin(-th)... position (L cos th, 0, L sin th) with th < 0 going down
        const Vec3 x = m.bodies[0].x;
        return std::atan2(x(2), x(0));
    };
    double c1, c2, c3, d1, d2, d3;
    const double t1 = pend(0.02, c1, d1), t2 = pend(0.01, c2, d2), t3 = pend(0.005, c3, d3);
    const double e1 = std::fabs(t1 - thref), e2 = std::fabs(t2 - thref), e3 = std::fabs(t3 - thref);
    check(e3 < 2.0e-4, fmt("angle at t = 1.2 s against the ODE: %.5f vs %.5f rad", t3, thref));
    check(std::log2(e1/e2) > 1.8 && std::log2(e2/e3) > 1.8, fmt("second order (orders %.2f %.2f)", std::log2(e1/e2), std::log2(e2/e3)));
    check(std::max(c1, std::max(c2, c3)) < 1.0e-10, fmt("constraint kept on position level: max |g| %.1e", std::max(c1, std::max(c2, c3))));
    check(d3 < 1.0e-3*mass*g*L, fmt("energy error %.1e J (of m g L = %.1f J)", d3, mass*g*L));
}

static void test_prismatic()
{
    std::cout<<"prismatic joints (float and spar)"<<std::endl;
    // float: prismatic to the ground along x (no rotation, no y/z motion); spar: prismatic to the
    // float along z, on a vertical spring to the ground; lateral force on the spar at an offset point
    model m;
    const double mA = 2.0, mB = 3.0, k = 30.0, F = 5.0;
    m.add_body(make_body(mA, Vec3(1,1,1), Vec3(0,0,0), "float"));
    m.add_body(make_body(mB, Vec3(2,2,1), Vec3(0,0,-2), "spar"));
    m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 0, Vec3(0,0,0), Vec3::UnitX()));
    joint &j = static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, 0, 1, Vec3(0,0,-1), Vec3::UnitZ())));
    Mat6 K = Mat6::Zero();
    K(2,2) = k;
    m.add_load(std::make_unique<linear6>(m, 1, Vec3(0,0,-2), K, Mat6::Zero()));
    m.add_load(std::make_unique<test_load>([&](const model &mm, double, VecX &f)
    {
        mm.add_force(f, 1, mm.bodies[1].point(Vec3(0.3,0,-1)), Vec3(F, 0, 0));
    }));
    m.bodies[1].x(2) = -2.0 + 0.2;    // axial offset 0.2: the spar oscillates along the joint
    timeint_settings s;
    timeint ti(s);
    const double tend = 1.0;
    run(m, ti, 0.0, tend, 0.002);
    const Vec3 xa = m.bodies[0].x, xb = m.bodies[1].x;
    const double xex = 0.5*F/(mA + mB)*tend*tend;
    check(std::fabs(xa(0) - xex) < 1.0e-9 && std::fabs(xb(0) - xex) < 1.0e-9,
          fmt("lateral: both bodies move with F/(mA+mB): x %.6f %.6f (exact %.6f)", xa(0), xb(0), xex));
    check((m.bodies[0].R - Mat3::Identity()).norm() < 1.0e-10 && (m.bodies[1].R - Mat3::Identity()).norm() < 1.0e-10, "no rotation");
    const double w = std::sqrt(k/mB);
    const double zex = 0.2*std::cos(w*tend);
    check(std::fabs(j.coordinate(m) - zex) < 1.0e-4, fmt("axial: joint coordinate %.5f vs %.5f (spar on its spring)", j.coordinate(m), zex));
    check(std::fabs(xa(2)) < 1.0e-10, fmt("float does not move axially (%.1e)", std::fabs(xa(2))));
    // the float is pushed by the spar through the joint: reaction on the float = mA a
    // force of the second joint only: -G_j^T lambda_j
    VecX g(5);
    MatX G = MatX::Zero(5, m.ndof());
    j.eval(m, g, G);
    const VecX Rj = -G.transpose()*ti.multipliers().tail(5);
    check(std::fabs(Rj(0) - mA*F/(mA + mB)) < 1.0e-8 && std::fabs(Rj(6) + mA*F/(mA + mB)) < 1.0e-8,
          fmt("joint force on the float %.6f N = mA F/(mA+mB) = %.6f N", Rj(0), mA*F/(mA + mB)));
}

static void test_fixed()
{
    std::cout<<"fixed joint"<<std::endl;
    model m;
    m.add_body(make_body(1.0, Vec3(0.1,0.1,0.1), Vec3(0,0,0)));
    m.add_body(make_body(2.0, Vec3(0.2,0.3,0.1), Vec3(1,0,0)));
    m.add_constraint(std::make_unique<joint>(m, joint_type::fixed, 0, 1, Vec3(0.5,0,0), Vec3::UnitZ()));
    m.add_load(std::make_unique<gravity_load>());
    m.bodies[0].v = Vec3(0, 0, 1.0);
    m.bodies[1].v = Vec3(0, 0, 1.0);
    const Mat3 Rrel0 = m.bodies[0].R.transpose()*m.bodies[1].R;
    timeint_settings s;
    timeint ti(s);
    run(m, ti, 0.0, 1.0, 0.01);
    const Vec3 cg = (m.bodies[0].x + 2.0*m.bodies[1].x)/3.0;
    const Vec3 cg_ex = Vec3(2.0/3.0, 0.0, 1.0 - 0.5*9.81);
    check((cg - cg_ex).norm() < 1.0e-9, fmt("common CoG falls freely (err %.1e)", (cg - cg_ex).norm()));
    check((m.bodies[0].R.transpose()*m.bodies[1].R - Rrel0).norm() < 1.0e-10 && std::fabs((m.bodies[1].x - m.bodies[0].x).norm() - 1.0) < 1.0e-10,
          "bodies stay rigidly connected");
    check(ti.stats.max_constraint < 1.0e-10, fmt("max |g| %.1e", ti.stats.max_constraint));
}

static void test_heavy_top()
{
    std::cout<<"heavy symmetric top on a spherical joint"<<std::endl;
    model m;
    const double mass = 1.0, l = 0.5, J3 = 0.02, w3 = 200.0;
    // axis tilted by 0.3 rad from the vertical (in the x-z plane), CoG at distance l from the pivot,
    // fast spin about the symmetry axis (J3 w3^2 >> 4 m g l J1)
    const double tilt = 0.3;
    rigid_body b = make_body(mass, Vec3(0.05, 0.05, J3), l*Vec3(std::sin(tilt), 0.0, std::cos(tilt)));
    b.R = Eigen::AngleAxisd(tilt, Vec3::UnitY()).toRotationMatrix();
    b.w = Vec3(0.0, 0.0, w3);
    m.add_body(b);
    m.add_load(std::make_unique<gravity_load>());
    m.add_constraint(std::make_unique<joint>(m, joint_type::spherical, -1, 0, Vec3::Zero(), Vec3::UnitZ()));
    timeint_settings s;
    s.rho_inf = 0.8;
    timeint ti(s);
    const double E0 = m.kinetic_energy() + m.potential_energy();
    // angular momentum about the pivot, vertical component (conserved)
    auto Lz = [&]() {const rigid_body &bb = m.bodies[0]; return (bb.angular_momentum() + bb.x.cross(bb.mass*bb.v))(2);};
    const double L0 = Lz();
    double prec = 0.0, last = 0.0;
    run(m, ti, 0.0, 2.0, 5.0e-4, [&](double)
    {
        const double a = std::atan2(m.bodies[0].x(1), m.bodies[0].x(0));
        double d = a - last;
        if(d > M_PI) d -= 2.0*M_PI;
        if(d < -M_PI) d += 2.0*M_PI;
        prec += d;
        last = a;
    });
    const double dE = std::fabs(m.kinetic_energy() + m.potential_energy() - E0)/E0;
    const double dL = std::fabs(Lz() - L0)/std::fabs(L0);
    check(dE < 1.0e-6, fmt("energy kept: rel. error %.1e", dE));
    check(dL < 2.0e-4, fmt("vertical angular momentum about the pivot kept: rel. error %.1e", dL));
    check(std::fabs(m.bodies[0].x.norm() - l) < 1.0e-10, fmt("pivot kept (|x| - l = %.1e)", m.bodies[0].x.norm() - l));
    // reference: Euler's equations about the fixed point, J_p dw + w x J_p w = l e3 x (R^T m g),
    // J_p = J + m l^2 (transverse), quaternion attitude, RK4 with small steps
    const Mat3 Jp = Vec3(0.05 + mass*l*l, 0.05 + mass*l*l, J3).asDiagonal();
    const Mat3 Jpi = Jp.inverse();
    Eigen::Quaterniond q(Eigen::AngleAxisd(tilt, Vec3::UnitY()));
    Vec3 w(0.0, 0.0, w3);
    auto rhs = [&](const Eigen::Quaterniond &qq, const Vec3 &ww, Eigen::Vector4d &dq, Vec3 &dw)
    {
        const Mat3 R = qq.normalized().toRotationMatrix();
        const Vec3 tau = Vec3(0,0,l).cross(R.transpose()*Vec3(0,0,-9.81*mass));
        dw = Jpi*(tau - ww.cross(Jp*ww));
        const Eigen::Quaterniond o(0.0, ww(0), ww(1), ww(2));
        const Eigen::Quaterniond d = qq*o;
        dq = 0.5*d.coeffs();
    };
    const double hr = 2.0e-6;
    for(int i=0; i<(int)std::lround(2.0/hr); ++i)
    {
        Eigen::Vector4d k1q, k2q, k3q, k4q; Vec3 k1w, k2w, k3w, k4w;
        auto add = [](const Eigen::Quaterniond &a_, const Eigen::Vector4d &d, double f){Eigen::Quaterniond r; r.coeffs() = a_.coeffs() + f*d; return r;};
        rhs(q, w, k1q, k1w);
        rhs(add(q, k1q, 0.5*hr), w + 0.5*hr*k1w, k2q, k2w);
        rhs(add(q, k2q, 0.5*hr), w + 0.5*hr*k2w, k3q, k3w);
        rhs(add(q, k3q, hr), w + hr*k3w, k4q, k4w);
        q.coeffs() += hr/6.0*(k1q + 2.0*k2q + 2.0*k3q + k4q);
        q.normalize();
        w += hr/6.0*(k1w + 2.0*k2w + 2.0*k3w + k4w);
    }
    const Vec3 xref = q.toRotationMatrix()*Vec3(0,0,l);
    const double ex = (m.bodies[0].x - xref).norm();
    check(ex < 2.0e-3*l, fmt("CoG at t = 2 s against Euler's equations about the pivot: error %.1e m (precession %.3f rad)", ex, prec));
}

static void test_cylindrical_planar()
{
    std::cout<<"cylindrical and planar joints"<<std::endl;
    {
        model m;
        rigid_body b = make_body(1.0, Vec3(0.1,0.1,0.2), Vec3(0,0,0));
        b.v = Vec3(0.3, 0.4, 1.0);
        b.w = Vec3(0.5, 0.2, 2.0);
        m.add_body(b);
        m.add_constraint(std::make_unique<joint>(m, joint_type::cylindrical, -1, 0, Vec3::Zero(), Vec3::UnitZ()));
        timeint ti;
        run(m, ti, 0.0, 1.0, 0.01);
        const rigid_body &r = m.bodies[0];
        check(std::fabs(r.x(0)) < 1.0e-10 && std::fabs(r.x(1)) < 1.0e-10 && std::fabs(r.x(2) - 1.0) < 1.0e-9,
              fmt("cylindrical: slides along z only (z = %.6f)", r.x(2)));
        check(std::fabs((r.R*Vec3::UnitZ()).dot(Vec3::UnitZ()) - 1.0) < 1.0e-12 && std::fabs(r.omega_global()(2) - 2.0) < 1.0e-9,
              fmt("cylindrical: rotates about z only (w_z = %.6f)", r.omega_global()(2)));
    }
    {
        model m;
        rigid_body b = make_body(1.0, Vec3(0.1,0.1,0.2), Vec3(0,0,0));
        b.v = Vec3(0.3, 0.4, 1.0);
        b.w = Vec3(0.5, 0.2, 2.0);
        m.add_body(b);
        m.add_constraint(std::make_unique<joint>(m, joint_type::planar, -1, 0, Vec3::Zero(), Vec3::UnitY()));
        m.add_load(std::make_unique<gravity_load>());
        timeint ti;
        run(m, ti, 0.0, 1.0, 0.01);
        const rigid_body &r = m.bodies[0];
        check(std::fabs(r.x(1)) < 1.0e-10 && std::fabs(r.x(0) - 0.3) < 1.0e-9 && std::fabs(r.x(2) - (1.0 - 0.5*9.81)) < 1.0e-9,
              fmt("planar (normal y): moves in the x-z plane (x %.6f, z %.6f)", r.x(0), r.x(2)));
        check(std::fabs(std::fabs((r.R*Vec3::UnitY()).dot(Vec3::UnitY())) - 1.0) < 1.0e-12, "planar: rotates about y only");
    }
}

static void test_tether()
{
    std::cout<<"tether (tension-only spring)"<<std::endl;
    // body hanging on a tether from the ground point (0,0,0), dropped from above the unstretched length
    model m;
    const double mass = 1.0, k = 400.0, L0 = 1.0;
    m.add_body(make_body(mass, Vec3(0.01,0.01,0.01), Vec3(0,0,-0.7)));
    m.add_load(std::make_unique<gravity_load>());
    spring_damper &sp = static_cast<spring_damper&>(m.add_load(std::make_unique<spring_damper>(m, -1, Vec3::Zero(), 0, Vec3(0,0,-0.7), k, 0.0, L0, 0.0, true)));
    m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 0, Vec3(0,0,-0.7), Vec3::UnitZ()));
    timeint_settings s;
    s.rho_inf = 1.0;
    timeint ti(s);
    const double E0 = m.kinetic_energy() + m.potential_energy();
    double Tmin = 1.0e30, Tmax = 0.0, zmin = 0.0;
    run(m, ti, 0.0, 3.0, 2.0e-4, [&](double){Tmin = std::min(Tmin, sp.tension(m)); Tmax = std::max(Tmax, sp.tension(m)); zmin = std::min(zmin, m.bodies[0].x(2));});
    const double dE = std::fabs(m.kinetic_energy() + m.potential_energy() - E0);
    // lowest point from energy: m g (0.7 - z) ... m g (h) = 1/2 k s^2 with h = drop to the lowest point
    // free fall 0.3 m, then the stretch s: m g (0.3 + s) = 1/2 k s^2
    const double g = 9.81, a = 0.5*k, bq = -mass*g, cq = -mass*g*0.3;
    const double sx = (-bq + std::sqrt(bq*bq - 4.0*a*cq))/(2.0*a);
    check(Tmin >= 0.0 && Tmax > 0.0, fmt("tension never negative (min %.2f N, max %.1f N)", Tmin, Tmax));
    check(std::fabs(-zmin - (L0 + sx)) < 2.0e-3, fmt("lowest point %.4f m vs energy balance %.4f m", -zmin, L0 + sx));
    check(dE < 1.0e-3*mass*g, fmt("energy through slack and taut phases: error %.1e J", dE));
}

static void test_newton_reuse()
{
    std::cout<<"modified Newton"<<std::endl;
    model m;
    rigid_body b = make_body(3.0, Vec3(1,1,1), Vec3(0.1,0,0));
    m.add_body(b);
    m.bodies[0].x = Vec3::Zero();
    Mat6 K = Mat6::Identity()*10.0, C = Mat6::Identity()*0.5;
    m.add_load(std::make_unique<linear6>(m, 0, Vec3::Zero(), K, C));
    m.bodies[0].x = Vec3(0.1,0,0);
    timeint ti;
    run(m, ti, 0.0, 10.0, 0.01);
    check(ti.stats.jacobians < ti.stats.steps/20, fmt("Jacobian rebuilt %.0f times in %.0f steps", (double)ti.stats.jacobians, (double)ti.stats.steps));
    check(ti.stats.failures==0, "no failed steps");
}

int main()
{
    std::cout<<coral::banner()<<std::endl;
    test_version();
    test_types();
    test_rotation();
    test_free_fall();
    test_oscillator();
    test_top();
    test_pendulum();
    test_prismatic();
    test_fixed();
    test_heavy_top();
    test_cylindrical_planar();
    test_tether();
    test_newton_reuse();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
