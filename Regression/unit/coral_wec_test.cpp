// Standalone tests of the REEF3D::CORAL wave-to-wire PTO (src/coral/pto): no MPI, no REEF3D.
// Build and run with coral_tests.sh (it builds the user controller coral_pto_user_example.c first and
// passes the library as the first argument):
//   ./coral_wec_test <path of coral_pto_user_example.so>
//
// Generator (coral_generator): ideal, pmsg and map against their formulas, gear ratio, motoring with
// reactive power; in a decay test the electrical energy is eta times the absorbed energy.
// Hydraulic PTO (coral_pto_hydraulic): steady piston speed against the analytical motor speed, motor
// pressure and valve losses; energy balance of the element (piston work = generator work + dissipation +
// gas and shaft energy); in a CORAL model under a sine force: the PTO's absorbed energy equals the
// element's piston work, the time step does not change the mean power, the checkpoint restores it.
// Link PTO (two points): along the line of a prismatic joint the same motion as the joint PTO; a flap on
// a hinge with an inclined actuator: q = change of the distance, the hinge moment F (r x e), absorbed
// energy = loss of mechanical energy.
// Controllers (coral_pto_wec_control): tuned passive and reactive gains of a 1-DOF absorber against the
// frequency domain optimum (reactive: F^2/(8 B)); the gain schedule switches the damping; the user library
// gives the same motion as the passive PTO with its gains and keeps its state over a checkpoint; MPC
// without constraints reaches the reactive optimum, with a force limit it keeps the limit; the real form
// of the radiation state space has the kernel of the complex one; the matrix exponential.
// Drive from a table (coral_joint_drive file): Hermite curve against the sampled sine, periodic, ramp; a
// body follows it. Keywords and channels of all of it.
// Architect: Hans Bihs
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/structure/coral_linear6.h"
#include"coral/hydro/coral_radiation.h"
#include"coral/hydro/coral_hydro_db.h"
#include"coral/hydro/coral_hydro_load.h"
#include"coral/pto/coral_pto.h"
#include"coral/pto/coral_pto_hydraulic.h"
#include"coral/pto/coral_generator.h"
#include"coral/pto/coral_pto_wec_control.h"
#include"coral/pto/coral_joint_drive.h"
#include"coral/pto/coral_joint_force.h"
#include"coral/io/coral_input.h"
#include"coral/core/coral_rotation.h"
#include<algorithm>
#include<fstream>
#include<iostream>
#include<string>
#include<cstdio>
#include<cmath>
#include<complex>
#include<memory>

using namespace coral;

static int nfail = 0;
static std::string user_lib;
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
template<class T> static T *find_load(model &m)
{
    for(auto &l : m.load_list)
        if(auto p = dynamic_cast<T*>(l.get())) return p;
    return nullptr;
}

// ---------------------------------------------------------------------------------------------------
static void test_generator()
{
    std::cout<<"generator"<<std::endl;
    generator_settings s;
    s.eta = 0.9; s.eta_conv = 0.95;
    generator_model g(s);
    const generator_power a = g.convert(-1000.0, 2.0), b = g.convert(1000.0, 2.0);
    check(std::fabs(a.P_grid - 2000.0*0.9*0.95) < 1.0e-9 && std::fabs(b.P_grid + 2000.0/0.9/0.95) < 1.0e-9,
          fmt("ideal: generating %.1f W, motoring %.1f W", a.P_grid, b.P_grid));
    // pmsg: losses 3/2 R i^2 + k_fe w^2, a held force costs power; gear ratio and efficiency
    s = generator_settings();
    s.kind = generator_kind::pmsg; s.kt = 50.0; s.R = 0.2; s.k_fe = 0.01; s.ratio = 10.0; s.eta_gear = 0.95;
    g = generator_model(s);
    const generator_power p = g.convert(-2000.0, 1.5);
    const double T = 2000.0/10.0*0.95, w = 15.0, I = T/50.0;
    const double Pe = T*w - 1.5*0.2*I*I - 0.01*w*w;
    const generator_power h = g.convert(-2000.0, 0.0);
    const generator_power mo = g.convert(2000.0, 1.5);
    check(std::fabs(p.w - w) < 1.0e-12 && std::fabs(p.T - T) < 1.0e-9 && std::fabs(p.I - I) < 1.0e-12 && std::fabs(p.P_el - Pe) < 1.0e-9
          && h.P_el < 0.0 && std::fabs(mo.P_el + 2000.0*1.5/0.95 + 1.5*0.2*std::pow(2000.0/10.0/0.95/50.0, 2) + 0.01*w*w) < 1.0e-9,
          fmt("pmsg with gear: shaft %.2f rad/s, %.1f N m, %.2f A, P_el %.2f W; holding costs, motoring pays the losses", p.w, p.T, p.I, p.P_el));
    // efficiency map, bilinear, end values outside
    std::ofstream("gen_map.dat")<<"# w \\ T\n0   100  200\n0   0.5  0.6\n10  0.7  0.9\n";
    s = generator_settings();
    std::string msg;
    const bool ok = generator_model::read_map("gen_map.dat", s, msg);
    std::remove("gen_map.dat");
    g = generator_model(s);
    const double e1 = g.map_eta(5.0, 150.0), e2 = g.map_eta(20.0, 300.0);
    check(ok && s.kind == generator_kind::map && std::fabs(e1 - 0.675) < 1.0e-12 && std::fabs(e2 - 0.9) < 1.0e-12,
          fmt("efficiency map: eta(5, 150) = %.4f, eta outside = %.3f", e1, e2));
    // in a decay test: electrical energy = eta eta_conv absorbed energy
    {
        coral_case c;
        c.read_string("time dt 0.001 end 1\ngravity off\nbody b mass 10 cog 0 0 0 offset 0.1 0 0\n"
                      "joint s prismatic ground b point 0 0 0 axis 1 0 0\nlinear6 k body b point 0 0 0 k 1000 0 0 0 0 0\n"
                      "pto p joint s c 20 generator ideal eta 0.9 eta_conv 0.95\n", ".");
        timeint ti;
        ti.initialize(c.m, 0.0);
        for(int i=0; i<2000; ++i)
            ti.step(c.m, i*1.0e-3, 1.0e-3);
        pto_load *pl = find_load<pto_load>(c.m);
        std::vector<std::string> ch;
        pl->channel_names(ch);
        check(std::fabs(pl->electrical_energy() - 0.9*0.95*pl->energy()) < 1.0e-9*pl->energy() && pl->energy() > 0.04
              && std::find(ch.begin(), ch.end(), "p_Pe") != ch.end(),
              fmt("decay: electrical energy %.6f J = eta eta_conv x absorbed %.6f J; channel p_Pe", pl->electrical_energy(), pl->energy()));
    }
}

// ---------------------------------------------------------------------------------------------------
static hydraulic_settings hyd_settings()
{
    hydraulic_settings s;
    s.A_a = s.A_b = 0.01;
    s.V_a0 = s.V_b0 = 0.01;
    s.A_v = 5.0e-4; s.dp_open = 2.0e4;
    s.V_hp = 0.1; s.p_hp_pre = 5.0e6; s.p_hp0 = 1.0e7;
    s.V_lp = 0.1; s.p_lp_pre = 5.0e5; s.p_lp0 = 1.0e6;
    s.D = 1.0e-4/(2.0*M_PI)*10.0;
    s.J = 1.0; s.b_gen = 5.0; s.eta_m = 0.9; s.c_w = 0.1;
    return s;
}

static void test_hydraulic()
{
    std::cout<<"hydraulic PTO"<<std::endl;
    // constant piston speed: Q = A v = D w, eta D dp = (b + c_w) w, valve loss rho/2 (Q/(Cd A_v))^2
    {
        hydraulic_settings s = hyd_settings();
        pto_hydraulic h(s);
        const double v = 0.5, dt = 0.01;
        pto_state st;
        st.qd = v;
        h.start(st);
        double F = 0.0;
        for(int i=1; i<=6000; ++i)
        {
            st.t = i*dt; st.q = v*st.t;
            F = h.force(st).F;
            h.accept(st, dt);
        }
        const auto y = h.accepted();
        const double ph = h.p_hp(y[2]), pl = h.p_lp(y[3]);
        const double Q = s.A_a*v, w = Q/s.D, dp = (s.b_gen + s.c_w)*w/(s.eta_m*s.D);
        const double dpv = 0.5*s.rho*std::pow(Q/(s.Cd*s.A_v), 2);
        const double Fe = -s.A_a*(ph + dpv) + s.A_a*(pl - dpv);
        check(std::fabs(y[4] - w) < 1.0e-6*w && std::fabs(ph - pl - dp) < 1.0e-6*dp,
              fmt("steady: motor %.4f rad/s (A v/D %.4f), dp %.4e Pa (%.4e)", y[4], w, ph - pl, dp));
        check(std::fabs(y[1] - ph - dpv) < 1.0e-3*dpv && std::fabs(pl - y[0] - dpv) < 1.0e-3*dpv && std::fabs(F - Fe) < 1.0e-5*std::fabs(Fe),
              fmt("valve loss %.1f Pa (rho/2 (Q/Cd A_v)^2 = %.1f), cylinder force %.1f N (%.1f)", y[1] - ph, dpv, F, Fe));
    }
    // energy balance of the element under a sine motion with leakage and the relief valve
    {
        hydraulic_settings s = hyd_settings();
        s.C_leak = 1.0e-11; s.p_relief = 1.2e7; s.K_r = 1.0e-9; s.kappa = 1.3;
        pto_hydraulic h(s);
        pto_state st;
        st.qd = 0.0;
        h.start(st);
        const double dt = 0.01, a = 0.4, w = 1.2;
        double Emax = 0.0;
        for(int i=1; i<=3000; ++i)
        {
            st.t = i*dt; st.q = a*std::sin(w*st.t); st.qd = a*w*std::cos(w*st.t);
            h.force(st);
            h.accept(st, dt);
            Emax = std::max(Emax, h.E_in);
        }
        const auto y = h.accepted();
        pto_hydraulic::state y0{};
        const double Es = h.gas_energy(y) + h.shaft_energy(y) - h.shaft_energy(y0);
        const double bal = h.E_in - h.E_gen - h.E_diss - Es;
        // the rest: the elastic energy of the oil in the chambers and the work p^2 dV/(2 beta) that the
        // constant bulk modulus leaves out of the balance (both about p/beta of the work)
        check(std::fabs(bal) < 3.0e-3*h.E_gen && h.E_diss > 0.0 && Emax > 0.0,
              fmt("energy balance: piston %.0f J = generator %.0f + dissipated %.0f + stored %.0f J", h.E_in, h.E_gen, h.E_diss, Es)
              + fmt(" (rest %.1e of the generator work)", bal/h.E_gen));
        // the sub-step: SDIRK2 at dt_sub 1e-3 s against 5e-5 s
        s.dt_sub = 5.0e-5;
        pto_hydraulic h2(s);
        st = pto_state();
        h2.start(st);
        for(int i=1; i<=3000; ++i)
        {
            st.t = i*dt; st.q = a*std::sin(w*st.t); st.qd = a*w*std::cos(w*st.t);
            h2.force(st);
            h2.accept(st, dt);
        }
        check(std::fabs(h2.E_in - h.E_in) < 1.0e-4*std::fabs(h.E_in) && std::fabs(h2.E_gen - h.E_gen) < 1.0e-4*h.E_gen,
              fmt("sub-step 1e-3 s against 5e-5 s: piston work %.2f J and %.2f J", h.E_in, h2.E_in));
    }
    // in a CORAL model: a float on a spring under a sine force, hydraulic PTO with a generator
    auto run = [](double dt, double tend, coral_case &c)
    {
        std::ofstream f("hyd_force.dat");
        for(int i=0; i<=400; ++i)
            f<<i*0.01*M_PI<<" "<<5.0e4*std::sin(2.0*i*0.01*M_PI)<<"\n";
        f.close();
        c.read_string("time dt 0.01 end 1\ngravity off\nbody b mass 10000 cog 0 0 0\njoint s prismatic ground b point 0 0 0 axis 0 0 1\n"
                      "linear6 k body b point 0 0 0 k 0 0 1e5 0 0 0\njforce f joint s file hyd_force.dat periodic ramp 3\n"
                      "pto p joint s hydraulic A_a 0.005 V_a0 0.005 A_v 3e-4 dp_open 2e4 V_hp 0.05 p_hp_pre 5e6 p_hp0 6e6 "
                      "V_lp 0.05 p_lp_pre 5e5 p_lp0 1e6 D 1.6e-5 J 1 b_gen 0.7 eta_m 0.9 c_w 0.01 "
                      "generator pmsg kt 2 R 0.05\n", ".");
        std::remove("hyd_force.dat");
        timeint ti;
        ti.initialize(c.m, 0.0);
        pto_load *pl = find_load<pto_load>(c.m);
        const int n = int(std::round(tend/dt)), n1 = int(std::round((tend - 4.0*M_PI)/dt));
        double E1 = 0.0;
        for(int i=0; i<n; ++i)
        {
            if(i == n1) E1 = pl->energy();
            ti.step(c.m, i*dt, dt);
        }
        return (pl->energy() - E1)/(4.0*M_PI);
    };
    {
        coral_case a, b;
        const double Pa = run(0.01, 30.0, a), Pb = run(0.005, 30.0, b);
        pto_load *pl = find_load<pto_load>(a.m);
        const pto_hydraulic *h = pl->hydraulic();
        check(h && std::fabs(pl->energy() - h->E_in) < 0.01*pl->energy() && pl->electrical_energy() > 0.0 && pl->electrical_energy() < h->E_gen,
              fmt("model: absorbed %.0f J = piston work %.0f J, electrical %.0f J < generator shaft %.0f J", pl->energy(), h->E_in,
                  pl->electrical_energy(), h->E_gen));
        check(Pa > 100.0 && std::fabs(Pa - Pb) < 0.01*Pa, fmt("mean absorbed power %.1f W at dt 0.01 s, %.1f W at 0.005 s", Pa, Pb));
        std::vector<std::string> ch;
        pl->channel_names(ch);
        bool chok = true;
        for(const char *x : {"p_pa", "p_pb", "p_php", "p_plp", "p_w", "p_Qm", "p_Pgen", "p_Pe", "p_Ee"})
            chok = chok && std::find(ch.begin(), ch.end(), x) != ch.end();
        std::vector<double> cv;
        pl->channel_values(a.m, 30.0, cv);
        check(chok && cv.size() == ch.size(), fmt("channels: %.0f names and values, pressures and the electrical power", ch.size()));
        // checkpoint: the whole model into a fresh one, the same force and energies
        coral_case c;
        run(0.01, 0.02, c);
        std::vector<double> v;
        a.m.save_state(v);
        size_t k = 0;
        const bool ok = c.m.load_state(v, k);
        pto_load *pc = find_load<pto_load>(c.m);
        check(ok && pc->force(c.m, 30.0) == pl->force(a.m, 30.0) && pc->hydraulic()->E_in == h->E_in
              && pc->electrical_energy() == pl->electrical_energy() && pc->hydraulic()->accepted() == h->accepted(),
              "checkpoint: pressures, accumulators, shaft and energies restored");
    }
}

// ---------------------------------------------------------------------------------------------------
static void test_link()
{
    std::cout<<"PTO between two points (link)"<<std::endl;
    // along the line of the prismatic joint: the same as the joint PTO
    {
        auto make = [](coral_case &c, bool link)
        {
            c.read_string(std::string("time dt 0.001 end 1\ngravity off\nbody b mass 10 cog 0 0 0 offset 0.1 0 0\n"
                          "joint s prismatic ground b point 0 0 0 axis 1 0 0\nlinear6 k body b point 0 0 0 k 1000 0 0 0 0 0\n")
                          + (link ? "pto p a ground -2 0 0 b b 0 0 0 c 20 coulomb 1\n" : "pto p joint s c 20 coulomb 1\n"), ".");
        };
        coral_case a, b;
        make(a, false);
        make(b, true);
        timeint ta, tb;
        ta.initialize(a.m, 0.0);
        tb.initialize(b.m, 0.0);
        double e = 0.0;
        for(int i=0; i<2000; ++i)
        {
            ta.step(a.m, i*1.0e-3, 1.0e-3);
            tb.step(b.m, i*1.0e-3, 1.0e-3);
            e = std::max(e, std::fabs(a.m.bodies[0].x(0) - b.m.bodies[0].x(0)));
        }
        check(e < 1.0e-10 && find_load<pto_load>(b.m)->is_link(), fmt("link along the joint axis = joint PTO: max difference %.1e m", e));
    }
    // a flap on a hinge (y axis at the origin), actuator from the ground (-1, 0, 0.5) to the flap (0, 0, 1)
    {
        coral_case c;
        c.read_string("time dt 0.001 end 1\ngravity off\nbody f mass 100 inertia 50 50 50 cog 0 0 1\n"
                      "joint h revolute ground f point 0 0 0 axis 0 1 0\nlinear6 k body f point 0 0 0 k 0 0 0 0 2000 0\n"
                      "pto p a ground -1 0 0.5 b f 0 0 1 c 300\n", ".");
        pto_load *pl = find_load<pto_load>(c.m);
        joint *jh = c.joints["h"];
        // rotate the flap by 0.2 rad about y: q = |p_b - p_a| - L0 and the hinge moment of a unit force
        const double th = 0.2;
        rigid_body &b = c.m.bodies[0];
        b.R = exp_so3(Vec3(0, th, 0));
        b.x = b.R*Vec3(0, 0, 1);
        const Vec3 pa(-1, 0, 0.5), pb = b.R*Vec3(0, 0, 1);
        const double q = (pb - pa).norm() - std::sqrt(1.25);
        VecX f = VecX::Zero(c.m.ndof());
        pl->apply(c.m, 1.0, f);
        const Vec3 M = b.R*f.segment<3>(3) + (b.x - Vec3::Zero()).cross(f.segment<3>(0));
        const Vec3 e = (pb - pa).normalized();
        const double My = pb.cross(e)(1);
        check(std::fabs(pl->coordinate(c.m) - q) < 1.0e-14 && std::fabs(M(1) - My) < 1.0e-12,
              fmt("flap at 0.2 rad: q %.6f m (|p_b - p_a| - L0 %.6f), hinge moment %.6f N m per N (r x e %.6f)", pl->coordinate(c.m), q, M(1), My));
        // decay: absorbed energy = loss of mechanical energy
        b.R = exp_so3(Vec3(0, 0.15, 0));
        b.x = b.R*Vec3(0, 0, 1);
        timeint ti;
        ti.initialize(c.m, 0.0);
        const double E0 = c.m.kinetic_energy() + c.m.potential_energy();
        for(int i=0; i<3000; ++i)
            ti.step(c.m, i*1.0e-3, 1.0e-3);
        const double E = c.m.kinetic_energy() + c.m.potential_energy();
        check(std::fabs((E0 - E) - pl->energy()) < 2.0e-3*E0 && pl->energy() > 0.3*E0 && jh,
              fmt("flap decay: absorbed %.4f J = loss of mechanical energy %.4f J", pl->energy(), E0 - E));
    }
}

// ---------------------------------------------------------------------------------------------------
// 1-DOF absorber: m q'' + B q' + C q = Fe sin(w t) + F_pto, mean absorbed power over 10 periods
static const double M1 = 100.0, C1 = 4000.0, B1 = 50.0, F1 = 100.0, W1 = 4.0;
static void write_excitation(const char *file)
{
    std::ofstream f(file);
    for(int i=0; i<=60000; ++i)
        f<<i*1.0e-3<<" "<<F1*std::sin(W1*i*1.0e-3)<<"\n";
}
static double absorber(const std::string &pto, double tend=40.0, double dt=0.005, coral_case *out=nullptr, double *Fmax=nullptr)
{
    coral_case own;
    coral_case &c = out ? *out : own;
    write_excitation("wec_exc.dat");
    std::ofstream("wec_tab.dat")<<"# w A B\n0.1 0 50\n20 0 50\n";
    c.read_string("time dt 0.01 end 1\ngravity off\nbody b mass 100 cog 0 0 0\njoint s prismatic ground b point 0 0 0 axis 1 0 0\n"
                  "linear6 k body b point 0 0 0 k 4000 0 0 0 0 0 c 50 0 0 0 0 0\njforce f joint s file wec_exc.dat\n" + pto + "\n", ".");
    timeint ti;
    ti.initialize(c.m, 0.0);
    pto_load *pl = find_load<pto_load>(c.m);
    const double Tp = 2.0*M_PI/W1, t1 = tend - 10.0*Tp;
    double E1 = 0.0, t = 0.0, F = 0.0;
    const int n = int(std::round(tend/dt));
    for(int i=0; i<n; ++i)
    {
        ti.step(c.m, i*dt, dt);
        t = (i + 1)*dt;
        if(E1 == 0.0 && t >= t1) E1 = pl->energy() - pl->power(c.m)*(t - t1);
        if(t >= t1) F = std::max(F, std::fabs(pl->force(c.m, t)));
    }
    if(Fmax) *Fmax = F;
    std::remove("wec_exc.dat");
    std::remove("wec_tab.dat");
    return (pl->energy() - E1)/(t - t1);
}

static std::string num(double x)
{
    char b[64];
    std::snprintf(b, sizeof(b), "%.17g", x);
    return b;
}

static void test_controllers()
{
    std::cout<<"WEC controllers"<<std::endl;
    const double Popt = F1*F1/(8.0*B1);
    // tuned reactive and passive
    {
        double k, c;
        pto_reactive(W1, M1, 0.0, B1, C1, 0.0, k, c);
        coral_case cc;
        const double Pr = absorber("pto p joint s tune reactive T " + num(2.0*M_PI/W1) + " table wec_tab.dat C 4000", 40.0, 0.005, &cc);
        const pto_load *pl = find_load<pto_load>(cc.m);
        check(std::fabs(k - (W1*W1*M1 - C1)) < 1.0e-9 && std::fabs(pl->settings.k - k) < 1.0e-6 && std::fabs(pl->settings.c - B1) < 1.0e-9
              && std::fabs(Pr - Popt) < 0.01*Popt,
              fmt("reactive: k %.1f N/m, c %.1f N s/m, mean power %.3f W (F^2/(8 B) = %.3f W)", pl->settings.k, pl->settings.c, Pr, Popt));
        const double Bopt = pto_bopt(W1, M1, 0.0, B1, C1), X = W1*M1 - C1/W1;
        const double v2 = F1*F1/((B1 + Bopt)*(B1 + Bopt) + X*X), Pp = 0.5*Bopt*v2;
        const double P = absorber("pto p joint s tune passive T " + num(2.0*M_PI/W1) + " table wec_tab.dat C 4000");
        check(std::fabs(P - Pp) < 0.01*Pp && P < Pr, fmt("passive: B_opt %.2f N s/m, mean power %.4f W (frequency domain %.4f W)", Bopt, P, Pp));
    }
    // gain schedule: damping 30 until 2 s, then 300 (linear in between from 2 to 2.5 s)
    {
        std::ofstream("wec_sched.dat")<<"0 0 30\n2 0 30\n2.5 0 300 5\n";
        coral_case c;
        absorber("pto p joint s schedule wec_sched.dat", 3.0, 0.005, &c);
        std::remove("wec_sched.dat");
        pto_load *pl = find_load<pto_load>(c.m);
        const double F = pl->force(c.m, 3.0), qd = pl->rate(c.m);
        check(std::fabs(F - (5.0 - 300.0*qd)) < 1.0e-9*std::max(1.0, std::fabs(F)), fmt("schedule: F %.4f N = F0 - c qd with c 300, F0 5 (%.4f)", F, 5.0 - 300.0*qd));
    }
    // user library: the same motion as the passive PTO with its gains; its state over a checkpoint
    if(user_lib.empty())
        check(false, "user library: no library given (run through coral_tests.sh)");
    else
    {
        coral_case a, b;
        const double Pa = absorber("pto p joint s k -2400 c 50", 20.0, 0.005, &a);
        const double Pb = absorber("pto p joint s user " + user_lib + " params 2 -2400 50", 20.0, 0.005, &b);
        check(std::fabs(a.m.bodies[0].x(0) - b.m.bodies[0].x(0)) < 1.0e-12 && Pa == Pb,
              fmt("user library: K, B from the library, the same motion (power %.6f W, %.6f W)", Pa, Pb));
        std::vector<double> v;
        pto_load *pb = find_load<pto_load>(b.m);
        pb->save_state(v);
        coral_case c;
        absorber("pto p joint s user " + user_lib + " params 2 -2400 50", 0.01, 0.005, &c);
        size_t k = 0;
        pto_load *pc = find_load<pto_load>(c.m);
        const bool ok = pc->load_state(v, k);
        std::vector<double> w;
        pc->save_state(w);
        check(ok && k == v.size() && v.back() == 4001.0 && w == v, fmt("user library: %.0f calls (start + every step) kept over the checkpoint", v.back()));
        bool thrown = false;
        try {coral_case d; d.read_string("time dt 0.01 end 1\nbody b mass 1\njoint s prismatic ground b point 0 0 0 axis 1 0 0\n"
                                         "pto p joint s user ./no_such_library.so\n", ".");}
        catch(const std::exception &e) {thrown = std::string(e.what()).find("cannot load") != std::string::npos;}
        check(thrown, "user library: a missing library stops the input with its name");
    }
    // MPC: without constraints the reactive optimum, with a force limit the limit kept
    {
        const double P = absorber("pto p joint s mpc Tc 0.05 N 60 C 4000 B 50 excitation wec_exc.dat");
        const double P2 = absorber("pto p joint s mpc Tc 0.05 N 60 C 4000 B 50 excitation wec_exc.dat", 40.0, 0.001);
        check(std::fabs(P - Popt) < 0.01*Popt && std::fabs(P2 - P) < 0.002*P,
              fmt("mpc without constraints: mean power %.3f W (optimum %.3f W), %.3f W at dt 0.001 s", P, Popt, P2));
        double Fm = 0.0;
        const double Pl = absorber("pto p joint s F_max 120 mpc Tc 0.05 N 60 C 4000 B 50 excitation wec_exc.dat", 40.0, 0.005, nullptr, &Fm);
        const double Pp = absorber("pto p joint s F_max 120 tune passive T " + num(2.0*M_PI/W1) + " table wec_tab.dat C 4000");
        check(Fm <= 120.0*(1.0 + 1.0e-9) && Pl < P && Pl > Pp, fmt("mpc with F_max 120 N: max |F| %.2f N, power %.3f W (passive with the limit %.3f W)", Fm, Pl, Pp));
        // stroke limit: the motion is kept near q_max
        coral_case c;
        double qmax = 0.0;
        {
            write_excitation("wec_exc.dat");
            c.read_string("time dt 0.01 end 1\ngravity off\nbody b mass 100 cog 0 0 0\njoint s prismatic ground b point 0 0 0 axis 1 0 0\n"
                          "linear6 k body b point 0 0 0 k 4000 0 0 0 0 0 c 50 0 0 0 0 0\njforce f joint s file wec_exc.dat\n"
                          "pto p joint s mpc Tc 0.05 N 60 C 4000 B 50 q_max 0.15 excitation wec_exc.dat\n", ".");
            std::remove("wec_exc.dat");
            timeint ti;
            ti.initialize(c.m, 0.0);
            for(int i=0; i<8000; ++i)
            {
                ti.step(c.m, i*0.005, 0.005);
                if(i > 4000) qmax = std::max(qmax, std::fabs(c.m.bodies[0].x(0)));
            }
        }
        // unconstrained the amplitude is F/(2 B w) = 0.25 m
        check(qmax < 0.16, fmt("mpc with q_max 0.15 m: max |q| %.4f m (0.25 m without)", qmax));
    }
    // real form of the radiation state space, matrix exponential
    {
        radiation_ss::pair q;
        q.lam.resize(2); q.b.resize(2); q.c.resize(2);
        q.lam << std::complex<double>(-0.5, 2.0), std::complex<double>(-1.5, 0.0);
        q.b << std::complex<double>(1.0, 0.3), std::complex<double>(0.7, 0.0);
        q.c << std::complex<double>(2.0, -1.0), std::complex<double>(-0.4, 0.0);
        MatX A;
        VecX B, C;
        radiation_ss::real_form(q, A, B, C);
        double e = 0.0;
        for(double t : {0.0, 0.3, 1.1, 2.5})
        {
            std::complex<double> k = 0.0;
            for(int r=0; r<2; ++r) k += q.c(r)*q.b(r)*std::exp(q.lam(r)*t);
            e = std::max(e, std::fabs(C.dot(coral_expm(A*t)*B) - k.real()));
        }
        check(e < 1.0e-12, fmt("radiation state space, real form: kernel error %.1e", e));
    }
}

// ---------------------------------------------------------------------------------------------------
// MPC of a heaving float with a hydrodynamic database: radiation memory as states (the database mode),
// the excitation predicted from CORAL's regular wave; the mean power against the bound |Fe|^2/(8 B(w))
static void test_mpc_hydro()
{
    std::cout<<"MPC with a hydrodynamic database"<<std::endl;
    const double rho = 1025.0, g = 9.81, Cz = rho*g*M_PI*4.0, w = 1.2, H = 1.0;
    auto Bw = [](double x) {const double r = x/1.5; return 2.0e4*r*r*std::exp(1.0 - r*r);};
    auto Xw = [&](double x) {return Cz*std::exp(-x*x/g*2.0);};
    hydro_db db;
    const int n = 80;
    db.resize(1, n, 1);
    db.rho = rho; db.g = g; db.depth = -1.0;
    db.body[0].name = "f";
    db.C(2,2) = Cz;
    db.Ainf(2,2) = 1.2e4;
    db.has_ainf = true;
    for(int k=0; k<n; ++k)
    {
        db.w[k] = 0.05*(k + 1);
        db.A[k](2,2) = 1.2e4 + 3.0e3*std::exp(-db.w[k]*db.w[k]);
        db.B[k](2,2) = Bw(db.w[k]);
        db.X[k][0](2) = Xw(db.w[k]);
    }
    db.write("wec_db.chd");
    auto run = [&](const std::string &pto, double &Fmax)
    {
        coral_case c;
        c.read_string("time dt 0.02 end 1\ngravity off\nbody f mass 25760 inertia 1e5 1e5 1e5 cog 0 0 0\n"
                      "joint s prismatic ground f point 0 0 0 axis 0 0 1\n"
                      "waves regular H 1 T " + num(2.0*M_PI/w) + "\nwaves ramp 10\n"
                      "hydro h file wec_db.chd bodies f state_space no_buoyancy\n" + pto + "\n", ".");
        timeint ti;
        ti.initialize(c.m, 0.0);
        pto_load *pl = find_load<pto_load>(c.m);
        const double dt = 0.02, tend = 120.0, t1 = tend - 5.0*2.0*M_PI/w;
        double E1 = 0.0;
        Fmax = 0.0;
        for(int i=0; i<int(std::round(tend/dt)); ++i)
        {
            ti.step(c.m, i*dt, dt);
            const double t = (i + 1)*dt;
            if(E1 == 0.0 && t >= t1) E1 = pl->energy() - pl->power(c.m)*(t - t1);
            if(t >= t1) Fmax = std::max(Fmax, std::fabs(pl->force(c.m, t)));
        }
        return (pl->energy() - E1)/(tend - t1);
    };
    const double Fe = 0.5*H*Xw(w), Pb = Fe*Fe/(8.0*Bw(w));
    double F1, F2;
    const double P = run("pto p joint s mpc Tc 0.1 N 100 hydro h mode 3", F1);
    const double Pl = run("pto p joint s F_max " + num(0.5*F1) + " mpc Tc 0.1 N 100 hydro h mode 3", F2);
    std::remove("wec_db.chd");
    check(P > 0.97*Pb && P < 1.01*Pb, fmt("heaving float, radiation memory as states: mean power %.1f W (bound |Fe|^2/(8 B) %.1f W)", P, Pb));
    check(F2 <= 0.5*F1*(1.0 + 1.0e-9) && Pl < P && Pl > 0.3*P, fmt("with half the force (%.0f N): %.1f W", 0.5*F1, Pl));
}

// ---------------------------------------------------------------------------------------------------
static void test_drive_table()
{
    std::cout<<"joint drive from a table"<<std::endl;
    model m;
    rigid_body rb;
    rb.mass = 2.0;
    rb.inertia = Vec3(1, 1, 1).asDiagonal();
    m.add_body(rb);
    m.gravity = Vec3::Zero();
    joint &j = static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, joint_type::prismatic, -1, 0, Vec3::Zero(), Vec3::UnitX())));
    const double a = 0.3, T = 2.0, w = 2.0*M_PI/T;
    drive_settings s;
    for(int i=0; i<=40; ++i)
    {
        s.t.push_back(i*0.05);
        s.q.push_back(a*std::sin(w*i*0.05));
    }
    s.periodic = true;
    s.ramp = 0.5;
    joint_drive &d = static_cast<joint_drive&>(m.add_load(std::make_unique<joint_drive>(j, false, s)));
    double eq = 0.0, ev = 0.0;
    for(double t=0.5; t<6.0; t+=0.0137)
    {
        double q, qd, qdd;
        d.prescribed(t, q, qd, qdd);
        eq = std::max(eq, std::fabs(q - a*std::sin(w*t)));
        ev = std::max(ev, std::fabs(qd - a*w*std::cos(w*t)));
    }
    // Catmull-Rom slopes: error about (w h)^2/6 of a w
    check(eq < 2.0e-4*a && ev < 6.0e-3*a*w, fmt("Hermite curve through the samples of a sine (40 per period, periodic): |q| error %.1e m, |qd| %.1e m/s", eq, ev));
    timeint ti;
    ti.initialize(m, 0.0);
    double e = 0.0;
    for(int i=0; i<2000; ++i)
    {
        ti.step(m, i*1.0e-3, 1.0e-3);
        double q, qd, qdd;
        d.prescribed((i + 1)*1.0e-3, q, qd, qdd);
        e = std::max(e, std::fabs(j.coordinate(m) - q));
    }
    check(e < 1.0e-4*a, fmt("the body follows the table: max |q - q_p| %.1e m", e));
    // ramp, scale, held ends, keyword
    std::ofstream("drive_tab.dat")<<"# t q\n0 0\n1 1\n2 0\n";
    coral_case c;
    c.read_string("time dt 0.01 end 1\nbody b mass 1 cog 0 0 0\njoint s prismatic ground b point 0 0 0 axis 0 0 1\n"
                  "drive d joint s file drive_tab.dat scale 2 ramp 1 q0 0.5\n", ".");
    std::remove("drive_tab.dat");
    const joint_drive *jd = find_load<joint_drive>(c.m);
    double q1, q2, q3, x;
    jd->prescribed(1.0, q1, x, x);
    jd->prescribed(3.0, q2, x, x);
    jd->prescribed(0.5, q3, x, x);
    double Q, Qd, Qdd;
    jd->table(0.5, Q, Qd, Qdd);
    check(jd && std::fabs(q1 - 2.5) < 1.0e-12 && std::fabs(q2 - 0.5) < 1.0e-12 && std::fabs(q3 - (0.5 + 2.0*0.5*Q)) < 1.0e-12,
          fmt("keyword drive file: q(1) %.3f, held after the table %.3f, ramp at 0.5 s %.4f", q1, q2, q3));
}

// ---------------------------------------------------------------------------------------------------
static void test_keywords()
{
    std::cout<<"keywords"<<std::endl;
    bool thrown = false;
    try {coral_case d; d.read_string("time dt 0.01 end 1\nbody b mass 1\npto p c 3\n", ".");}
    catch(const std::exception &e) {thrown = std::string(e.what()).find("'joint <joint>' or 'a <body>") != std::string::npos;}
    check(thrown, "pto without joint or points: the message names both");
    thrown = false;
    try {coral_case d; d.read_string("time dt 0.01 end 1\nbody b mass 1\njoint s prismatic ground b point 0 0 0 axis 1 0 0\n"
                                     "pto p joint s hydraulic A_a 0.01 V_a0 0.01\n", ".");}
    catch(const std::exception &e) {thrown = std::string(e.what()).find("'A_v' is required") != std::string::npos;}
    check(thrown, "hydraulic: a missing parameter is named");
    thrown = false;
    try {coral_case d; d.read_string("time dt 0.01 end 1\nbody b mass 1\njoint s prismatic ground b point 0 0 0 axis 1 0 0\n"
                                     "pto p joint s latching 1 1e5 mpc Tc 0.1 excitation x.dat\n", ".");}
    catch(const std::exception &e) {thrown = std::string(e.what()).find("not with latching") != std::string::npos;}
    check(thrown, "mpc with latching is refused");
    coral_case c;
    c.read_string("time dt 0.01 end 1\nbody b mass 1 cog 0 0 0\njoint s prismatic ground b point 0 0 0 axis 1 0 0\n"
                  "pto p a ground 0 0 -3 b b 0 0 0 c 10 generator pmsg kt 20 R 0.1 k_fe 0.001 ratio 5 eta_gear 0.97 eta_conv 0.96\n", ".");
    const pto_load *pl = find_load<pto_load>(c.m);
    std::vector<std::string> ch;
    pl->channel_names(ch);
    check(pl->is_link() && pl->has_generator() && pl->settings.generator.kind == generator_kind::pmsg && pl->settings.generator.ratio == 5.0
          && pl->settings.generator.eta_conv == 0.96 && ch.back() == "p_I",
          "pto link with a pmsg generator, gear and converter: settings and channels");
}

int main(int argc, char **argv)
{
    if(argc > 1) user_lib = argv[1];
    test_generator();
    test_hydraulic();
    test_link();
    test_controllers();
    test_mpc_hydro();
    test_drive_table();
    test_keywords();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
