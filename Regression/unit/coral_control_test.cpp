// REEF3D::CORAL unit test of the wind turbine controller (coral/control/coral_turbine_controller) and the rotor
// tuning (coral/aero/coral_rotor_tune): the torque law of the NREL 5 MW baseline continuous through its regions, the
// speed filter's step response, the PI pitch loop with gain scheduling, rate limit and anti-windup, the PI gains of the
// NREL 5 MW from its published sensitivity and inertia; the tuning of a Glauert-optimum rotor (Cp_max, K, rated
// power held above rated); the closed loop in the multibody model: tip speed ratio lambda_opt below rated, rated speed,
// power and the tuned pitch above it; the keyword and its channels; the controller state round trip; the floating
// compensation: the platform's part of the speed out of the pitch loop, the sensitivities, a rotor on a pitching platform.
// Architect: Hans Bihs
#include"coral/control/coral_turbine_controller.h"
#include"coral/aero/coral_rotor.h"
#include"coral/aero/coral_rotor_tune.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/io/coral_input.h"
#include"coral/wind/coral_wind.h"
#include<cmath>
#include<cstdio>
#include<fstream>
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
static const double deg = M_PI/180.0;
static const double U0 = 10.0, R = 40.0, LAM = 7.0, OM = LAM*U0/R;
static airfoil thin()
{
    airfoil af;
    af.tables.emplace_back();
    polar_table &p = af.tables.back();
    for(int i=-180; i<=180; ++i)
    {
        const double a = i*deg;
        // 2 pi alpha up to 12 deg, a flat plate beyond (smooth enough for the tests)
        const double cl = std::fabs(a) < 12*deg ? 2*M_PI*a : 2*M_PI*12*deg*(a > 0 ? 1 : -1)*std::cos(a)*std::cos(a)/std::pow(std::cos(12*deg), 2);
        p.alpha.push_back(a); p.cl.push_back(cl); p.cd.push_back(0.008 + 1.2*std::sin(a)*std::sin(a)); p.cm.push_back(-0.05*std::sin(a));
    }
    return af;
}

static rotor_spec design(int n=30)
{
    rotor_spec sp;
    sp.B = 3;
    sp.hub = Vec3(0, 0, 90.0);
    sp.axis = Vec3(1, 0, 0);
    sp.R_hub = 2.0;
    sp.airfoils.push_back(thin());
    const double ad = 6.0*deg, cld = 2*M_PI*ad;
    for(int i=0; i<=n; ++i)
    {
        const double r = sp.R_hub + (R - sp.R_hub)*i/n;
        const double phi = 2.0/3.0*std::atan(1.0/(LAM*r/R));
        blade_station s;
        s.r = r;
        s.chord = std::min(8.0*M_PI*r*(1.0 - std::cos(phi))/(sp.B*cld), 4.0);
        s.twist = phi - ad;
        sp.stations.push_back(s);
    }
    sp.hub_mass = 50000.0;
    sp.hub_J = 1.0e5;
    sp.structure.frac = {0.0, 0.5, 1.0};
    sp.structure.mass = {400.0, 250.0, 100.0};
    sp.structure.pitch_axis = sp.structure.twist = sp.structure.EI_flap = sp.structure.EI_edge = {0, 0, 0};
    sp.omega0 = OM;
    sp.aero.rho = 1.225;
    return sp;
}

static double sp_twist(const rotor_spec &sp, int i) {return sp.stations[i].twist;}

static rotor_def steady_def(const rotor_spec &sp)
{
    rotor_def rd;
    rd.B = sp.B; rd.R_hub = sp.R_hub; rd.R_tip = sp.stations.back().r;
    rd.stations = sp.stations; rd.airfoils = sp.airfoils;
    return rd;
}

static void test_torque_law()
{
    std::cout<<"torque law (NREL 5 MW baseline)"<<std::endl;
    turbine_controller c;
    const turbine_controller_settings &s = c.set;
    const double e = 1e-9;
    auto T = [&](double w) {return c.torque_law(w, 0.0);};
    const double wt = s.w_trans();
    check(T(s.w_cutin)==0.0 && std::fabs(T(s.w_rgn2 - e) - T(s.w_rgn2 + e)) < 1e-3 && std::fabs(T(wt - e) - T(wt + e)) < 1e-3
          && std::fabs(T(s.w_rated - e) - T(s.w_rated + e)) < 1e-3,
          fmt("continuous at cut-in, region 2 (%.1f rpm), the transition (%.1f rpm) and rated (%.1f rpm)", s.w_rgn2*30/M_PI, wt*30/M_PI, s.w_rated*30/M_PI));
    check(std::fabs(T(100.0) - s.K*1e4) < 1e-9 && std::fabs(T(130.0) - s.P_rated/130.0) < 1e-9 && std::fabs(c.torque_law(100.0, 0.02) - s.P_rated/100.0) < 1e-9,
          "region 2: K w^2; region 3: P_rated/w, also below rated speed with pitch above pitch_rgn3");
    check(std::fabs(T(s.w_rated) - s.P_rated/s.w_rated) < 1e-6 && T(s.w_rated) < s.T_max,
          fmt("rated torque %.0f N m below T_max %.0f N m", T(s.w_rated), s.T_max));
}

static void test_filter_pi()
{
    std::cout<<"speed filter, PI pitch"<<std::endl;
    turbine_controller_settings s;
    s.torque_control = false;
    turbine_controller c(s);
    c.initialize(0.0, 100.0, 0.0, 0.0);
    const double dt = 0.01;
    for(int i=1; i<=100; ++i) c.update(i*dt, 120.0);
    const double ex = 120.0 - 20.0*std::exp(-s.omega_c*1.0);
    check(std::fabs(c.filtered_speed() - ex) < 1e-9, fmt("filter: 100 -> 120 rad/s step, %.6f after 1 s = 120 - 20 exp(-omega_c t)", c.filtered_speed()));
    // PI with the filter fast, no gain scheduling: theta = KP e + KI e t for a constant error
    s.omega_c = 1.0e9; s.theta_K = 1.0e9; s.theta_rate = 10.0;
    turbine_controller p(s);
    p.initialize(0.0, s.w_ref, 0.0, 0.0);
    const double err = 0.5;
    for(int i=1; i<=200; ++i) p.update(i*dt, s.w_ref + err);
    check(std::fabs(p.pitch() - (s.KP*err + s.KI*err*2.0)) < 1e-12, fmt("PI: pitch %.6f rad = KP e + KI e t after 2 s", p.pitch()));
    // gain scheduling: GK = 1/(1 + theta/theta_K) with the last command
    turbine_controller_settings g = s;
    g.theta_K = 0.1;
    turbine_controller q(g);
    q.initialize(0.0, g.w_ref, 0.1, 0.0);
    q.update(dt, g.w_ref + err);
    const double GK = 0.5, I = 0.1/(GK*g.KI) + err*dt;
    check(std::fabs(q.pitch() - GK*(g.KP*err + g.KI*I)) < 1e-12, "gain scheduling GK = 1/(1 + theta/theta_K), integrator started at the measured pitch");
    // rate limit and anti-windup
    turbine_controller r(s);
    r.initialize(0.0, s.w_ref, 0.0, 0.0);
    r.update(dt, s.w_ref + 1000.0);
    check(std::fabs(r.pitch() - s.theta_rate*dt) < 1e-12, "rate limit");
    turbine_controller w(s);
    w.initialize(0.0, s.w_ref, 0.0, 0.0);
    for(int i=1; i<=500; ++i) w.update(i*dt, s.w_ref - 10.0);       // below the reference for 5 s: held at 0
    w.update(501*dt, s.w_ref + err);
    check(w.pitch() > 0.0 && w.integral() <= err*dt + 1e-15, fmt("anti-windup: after 5 s below the reference the pitch moves at once (%.2e rad)", w.pitch()));
    std::vector<double> st;
    p.save_state(st);
    turbine_controller p2(s);
    size_t k = 0;
    check(p2.load_state(st, k) && k==st.size() && p2.pitch()==p.pitch() && p2.integral()==p.integral() && p2.filtered_speed()==p.filtered_speed(),
          "state round trip (checkpoint)");
}

static void test_gains()
{
    std::cout<<"PI gains"<<std::endl;
    // NREL 5 MW: rotor 38 677 040 kg m^2 + 97^2 534.116, 12.1 rpm, dP/dtheta(0) = -25.52 MW/rad, 0.6 rad/s, 0.7
    double KP, KI;
    turbine_controller_settings::pi_gains(38677040.0 + 97.0*97.0*534.116, 12.1*M_PI/30.0, 97.0, -25.52e6, 0.6, 0.7, KP, KI);
    check(std::fabs(KP/0.01882681 - 1.0) < 3e-3 && std::fabs(KI/0.008068634 - 1.0) < 3e-3,
          fmt("NREL 5 MW: KP %.6f s (published 0.018827), KI %.6f (0.008069)", KP, KI));
}

static void test_tuning()
{
    std::cout<<"rotor tuning"<<std::endl;
    const rotor_spec sp = design();
    const rotor_def rd = steady_def(sp);
    const double Pr = 3.0e6, Or = 2.0;
    const rotor_tuning tu = tune_rotor(rd, 1.225, Pr, Or);
    const double cpl = rotor_steady(rd, 10.0, (tu.lambda_opt - 0.2)*10.0/R, 0.0, 1.225).Cp, cph = rotor_steady(rd, 10.0, (tu.lambda_opt + 0.2)*10.0/R, 0.0, 1.225).Cp;
    check(cpl < tu.cp_max && cph < tu.cp_max && std::fabs(tu.K - 0.5*1.225*M_PI*std::pow(R, 5)*tu.cp_max/std::pow(tu.lambda_opt, 3)) < 1e-6*tu.K,
          fmt("Cp_max %.4f at tip speed ratio %.3f, K = rho pi R^5 Cp/(2 lambda^3) = %.0f N m s^2", tu.cp_max, tu.lambda_opt, tu.K));
    double e = 0.0, thp = -1.0;
    bool mono = true;
    for(const auto &p : tu.curve)
        if(p.U > tu.U_rated)
        {
            e = std::max(e, std::fabs(p.P/Pr - 1.0));
            mono = mono && p.theta > thp; thp = p.theta;
        }
    check(e < 1e-6 && mono && tu.dPdtheta0 < 0.0 && tu.theta_K > 0.0,
          fmt("rated wind %.2f m/s; above it the pitch holds 3 MW (%.1e), rising with the wind; dP/dtheta(0) %.2f MW/rad, theta_K %.2f deg",
              tu.U_rated, e, tu.dPdtheta0/1e6, tu.theta_K/deg));
}

// a rotor with the controller in the model, run to a steady state in steady wind
struct closed_loop
{
    model m;
    rotor_handles h;
    controller_load *ctrl = nullptr;
    rotor_tuning tu;
    closed_loop(double U, double rpm0)
    {
        rotor_spec sp = design(15);
        sp.omega0 = rpm0*M_PI/30.0;
        sp.dt.N = 50.0; sp.J_gen = 20.0;
        sp.dt.law = generator_law::command;
        m.gravity = Vec3::Zero();
        h = build_rotor(m, std::make_shared<wind_field>(U, 90.0), sp);
        tu = tune_rotor(steady_def(sp), 1.225, 3.0e6, 2.0);
        turbine_controller_settings cs;
        const double N = 50.0;
        cs.w_ref = N*2.0; cs.w_rated = 0.99*cs.w_ref; cs.w_cutin = 0.571*cs.w_ref; cs.w_rgn2 = 0.742*cs.w_ref;
        cs.P_rated = 3.0e6; cs.K = tu.K/(N*N*N); cs.T_max = 1.1*cs.P_rated/cs.w_rated; cs.T_rate = 0.35*cs.P_rated/cs.w_rated;
        cs.theta_K = tu.theta_K;
        turbine_controller_settings::pi_gains(h.J + 20.0*N*N, 2.0, N, tu.dPdtheta0, 0.6, 0.7, cs.KP, cs.KI);
        controller_load::io f;
        drivetrain *dt = h.drive;
        rotor_aero *ra = h.aero;
        f.speed = [dt](const model &mm) {return dt->omega_hss(mm);};
        f.torque = [dt](const model &mm) {return dt->generator_torque(mm);};
        f.set_torque = [dt](double T) {dt->set_torque_command(T);};
        f.pitch = [ra](const model&) {return ra->pitch(0, 0.0);};
        f.set_pitch = [ra](double th) {for(int q=0; q<3; ++q) ra->set_pitch(q, th);};
        ctrl = &static_cast<controller_load&>(m.add_load(std::make_unique<controller_load>(cs, f, "ctrl")));
    }
    void run(double T)
    {
        timeint ti;
        ti.initialize(m, 0.0);
        const double dt = 0.02;
        for(int i=0; i<(int)(T/dt); ++i) ti.step(m, i*dt, dt);
    }
};

static void test_closed_loop()
{
    std::cout<<"closed loop"<<std::endl;
    closed_loop a(10.0, 14.0);
    a.run(200.0);
    const double lam = a.h.drive->omega_rotor(a.m)*R/10.0;
    check(std::fabs(lam/a.tu.lambda_opt - 1.0) < 2e-3 && a.ctrl->controller().region()==2 && a.ctrl->controller().pitch()==0.0,
          fmt("10 m/s: region 2, tip speed ratio %.4f = lambda_opt %.4f, pitch 0", lam, a.tu.lambda_opt));
    closed_loop b(17.0, 18.0);
    b.run(200.0);
    double th15 = 0.0;
    for(const auto &p : b.tu.curve) if(std::fabs(p.U - 17.0) < 1e-9) th15 = p.theta;
    std::vector<station_load> L;
    b.h.aero->station_loads(b.m, 200.0, L);
    const auto v = b.h.aero->values(b.m, 200.0, L);
    check(std::fabs(b.h.drive->omega_rotor(b.m) - 2.0) < 1e-4 && std::fabs(b.ctrl->controller().pitch() - th15) < 1e-4 && std::fabs(v.P/3.0e6 - 1.0) < 1e-3
          && b.ctrl->controller().region()==3,
          fmt("17 m/s: region 3, %.5f rad/s (2), pitch %.4f deg = the tuning's %.4f deg, power %.1f kW", b.h.drive->omega_rotor(b.m),
              b.ctrl->controller().pitch()/deg, th15/deg, v.P/1e3));
}

static void test_input()
{
    std::cout<<"keyword"<<std::endl;
    {
        std::ofstream f("ctrl_test_af.dat");
        const airfoil af = thin();
        for(size_t i=0; i<af.tables[0].alpha.size(); ++i)
            f<<af.tables[0].alpha[i]/deg<<" "<<af.tables[0].cl[i]<<" "<<af.tables[0].cd[i]<<"\n";
        std::ofstream b("ctrl_test_blade.dat");
        const rotor_spec sp = design(10);
        b<<"test\n"<<sp.stations.size()<<"   NumBlNds\nBlSpn BlCrvAC BlSwpAC BlCrvAng BlTwist BlChord BlAFID\n";
        for(const blade_station &s : sp.stations)
            b<<s.r - 2.0<<" 0 0 0 "<<s.twist/deg<<" "<<s.chord<<" 1\n";
    }
    coral_case c;
    c.read_string("time dt 0.02 end 1\n"
                  "wind steady U 15 z_ref 90\n"
                  "rotor rot hub 0 0 90 axis 1 0 0 R_hub 2 blade ctrl_test_blade.dat airfoils 1 ctrl_test_af.dat hub_mass 50000 hub_inertia 1e7 \\\n"
                  "   rpm 19 gearbox 50 gen_inertia 20 pitch_joints\n"
                  "controller ctrl rotor rot tune P_rated 3e6 rpm_rated 19.0986 write ctrl_test_tuning.dat KP 0.05\n", ".");
    std::ifstream tf("ctrl_test_tuning.dat");
    const bool wrote = (bool)tf;
    std::remove("ctrl_test_af.dat"); std::remove("ctrl_test_blade.dat"); std::remove("ctrl_test_tuning.dat");
    controller_load *cl = nullptr;
    for(auto &l : c.m.load_list)
        if(l->name()=="ctrl") cl = dynamic_cast<controller_load*>(l.get());
    check(cl && wrote && cl->controller().set.KP==0.05 && std::fabs(cl->controller().set.w_ref - 50.0*2.0) < 1e-3
          && c.rotors.at("rot").drive->settings.law==generator_law::command,
          "controller with tune (operating curve written), KP given over the tuning, the generator on command");
    std::vector<std::string> n;
    if(cl) cl->channel_names(n);
    check(n.size()==4 && n[0]=="ctrl_wf" && n[2]=="ctrl_pitch", "channels ctrl_wf _T _pitch _region");
    auto err = [](const std::string &s)
    {
        try {coral_case k; k.read_string(s, ".");} catch(const std::exception &e) {return std::string(e.what());}
        return std::string();
    };
    check(err("controller c rotor nope\n").find("unknown rotor") != std::string::npos, "unknown rotor: error");
}

// the floating compensation: the controller alone, the sensitivities, and a rotor on a pitching platform
static void write_rotor_files(int n)
{
    std::ofstream f("ctrl_test_af.dat");
    const airfoil af = thin();
    for(size_t i=0; i<af.tables[0].alpha.size(); ++i)
        f<<af.tables[0].alpha[i]/deg<<" "<<af.tables[0].cl[i]<<" "<<af.tables[0].cd[i]<<"\n";
    std::ofstream b("ctrl_test_blade.dat");
    const rotor_spec sp = design(n);
    b<<"test\n"<<sp.stations.size()<<"   NumBlNds\nBlSpn BlCrvAC BlSwpAC BlCrvAng BlTwist BlChord BlAFID\n";
    for(const blade_station &st : sp.stations)
        b<<st.r - 2.0<<" 0 0 0 "<<st.twist/deg<<" "<<st.chord<<" 1\n";
    // blade mass (9.5 t per blade as validations 19, 21 - 23): the rotor's inertia sets the pitch loop's dynamics
    std::ofstream e("ctrl_test_ed.dat");
    e<<"---- ElastoDyn blade ----\ntest\n---- BLADE PARAMETERS ----\n3   NBlInpSt\n1.0   AdjBlMs\n1.0   AdjFlSt\n1.0   AdjEdSt\n"
       "---- DISTRIBUTED BLADE PROPERTIES ----\nBlFract PitchAxis StrcTwst BMassDen FlpStff EdgStff\n(-) (-) (deg) (kg/m) (Nm^2) (Nm^2)\n"
       "0.0 0.25 0 400 1e10 1e10\n0.5 0.25 0 250 1e9 2e9\n1.0 0.25 0 100 1e7 2e7\n";
}

static void test_floating()
{
    std::cout<<"floating compensation"<<std::endl;
    // the controller alone in region 3: the generator speed that the drivetrain model predicts for the nacelle's motion
    // (the exact solution of dw/dt = A w + B v for v = v0 sin(W t)) leaves the pitch alone; without the compensation the
    // pitch answers it
    {
        turbine_controller_settings cs;
        cs.float_A = -0.03; cs.float_B = -2.0;
        const double v0 = 0.5, W = 0.33, A = cs.float_A, Bf = cs.float_B;
        auto wp = [&](double t)
        {
            // particular solution plus the transient that starts from 0
            const double d = A*A + W*W;
            const double ps = Bf*v0*(-A*std::sin(W*t) - W*std::cos(W*t))/d;
            return ps + Bf*v0*W/d*std::exp(A*t);
        };
        double r[2];
        for(int k=0; k<2; ++k)
        {
            turbine_controller c(cs);
            if(k==1) c.set.float_B = 0.0;
            c.initialize(0.0, cs.w_ref, 0.15, 0.0);
            double lo = 1e9, hi = -1e9;
            for(int i=1; i<=6000; ++i)
            {
                const double t = 0.01*i;
                c.update(t, cs.w_ref + wp(t), v0*std::sin(W*(t - 0.005)));
                lo = std::min(lo, c.pitch()); hi = std::max(hi, c.pitch());
            }
            r[k] = hi - lo;
        }
        check(r[0] < 2e-3*r[1], fmt("controller: the platform's part of the speed out of the pitch loop: pitch range %.2e rad "
                                    "with the compensation, %.3f rad without", r[0], r[1]));
    }
    // the sensitivities above rated: thrust rises with the wind at a fixed pitch but falls at a fixed torque (the
    // cause of the negative damping); the pitch on the operating curve
    const rotor_spec sp = design();
    const rotor_def rd = steady_def(sp);
    const rotor_tuning tu = tune_rotor(rd, 1.225, 3.0e6, 2.0);
    const rotor_sensitivity rs = rotor_sensitivities(rd, tu, 1.225, 3.0e6, 2.0, 17.0);
    double th17 = 0.0;
    for(const auto &p : tu.curve) if(std::fabs(p.U - 17.0) < 1e-9) th17 = p.theta;
    const double dTdU_Q = rs.dT_dU - rs.dT_dtheta*rs.dQ_dU/rs.dQ_dtheta;
    const double d = 1e-4;
    const double fd = (rotor_steady(rd, 17.0, 2.0, rs.theta + d, 1.225).Q - rotor_steady(rd, 17.0, 2.0, rs.theta - d, 1.225).Q)/(2*d);
    check(std::fabs(rs.theta - th17) < 1e-8 && std::fabs(rs.Q*2.0/3.0e6 - 1.0) < 1e-6 && rs.dT_dU > 0.0 && dTdU_Q < 0.0
          && std::fabs(rs.dQ_dtheta/fd - 1.0) < 1e-3,
          fmt("17 m/s: pitch %.3f deg of the curve, dT/dU %.0f N/(m/s) at a fixed pitch, %.0f at a fixed torque, dQ/dtheta %.3e",
              rs.theta/deg, rs.dT_dU, dTdU_Q, rs.dQ_dtheta));
    // the rotor on a platform pitching about a hinge 90 m below the hub (19 s), steady wind of 17 m/s, the onshore
    // gains (0.6 rad/s): the thrust step starts the platform pitching; without the compensation the oscillation grows,
    // with it it decays
    write_rotor_files(12);
    double amp[2][2];
    bool kw = false;
    for(int k=0; k<2; ++k)
    {
        coral_case c;
        c.read_string(std::string("time dt 0.05 end 1\ngravity off\nwind steady U 17 z_ref 90\n"
                      "body plat mass 1e5 cog 0 0 90 inertia 1e8 3.2e9 1e8\n"
                      "joint hinge revolute ground plat point 0 0 0 axis 0 1 0\n"
                      "linear6 k body plat point 0 0 0 k 0 0 0 0 4.4e8 0 c 0 0 0 0 2.4e7 0\n"
                      "rotor rot hub 0 0 90 axis 1 0 0 R_hub 2 blade ctrl_test_blade.dat airfoils 1 ctrl_test_af.dat structure ctrl_test_ed.dat "
                      "hub_mass 50000 hub_inertia 1e5 nacelle plat rpm 19.0986 pitch 10 gearbox 50 gen_inertia 20\n"
                      "controller ctrl rotor rot tune P_rated 3e6 rpm_rated 19.0986") + (k ? " floating\n" : "\n"), ".");
        const int b = c.m.find_body("plat");
        if(k==1)
            for(auto &l : c.m.load_list)
                if(auto *cl = dynamic_cast<controller_load*>(l.get()))
                {
                    std::vector<std::string> n;
                    cl->channel_names(n);
                    const auto &st = cl->controller().set;
                    kw = n.size()==6 && n[4]=="ctrl_wplat" && n[5]=="ctrl_vnac" && st.float_A < 0.0 && st.float_B < 0.0 && st.float_gain==1.0;
                }
        timeint ti;
        ti.initialize(c.m, 0.0);
        double lo[2] = {1e9, 1e9}, hi[2] = {-1e9, -1e9};
        for(int i=0; i<6000; ++i)
        {
            ti.step(c.m, i*0.05, 0.05);
            const double t = (i + 1)*0.05;
            const double ang = std::atan2(-(c.m.body(b).R(2, 0)), c.m.body(b).R(0, 0));
            const int w = t > 40.0 && t <= 100.0 ? 0 : t > 240.0 ? 1 : -1;
            if(w >= 0) {lo[w] = std::min(lo[w], ang); hi[w] = std::max(hi[w], ang);}
        }
        amp[k][0] = 0.5*(hi[0] - lo[0])/deg; amp[k][1] = 0.5*(hi[1] - lo[1])/deg;
    }
    std::remove("ctrl_test_af.dat"); std::remove("ctrl_test_blade.dat"); std::remove("ctrl_test_ed.dat");
    check(kw, "keyword controller ... tune floating: float_A, float_B from the BEM (both negative), channels _wplat _vnac");
    check(amp[0][1] > amp[0][0] && amp[1][1] < 0.2*amp[1][0],
          fmt("rotor on a pitching platform: amplitude %.2f -> %.2f deg without the compensation, %.2f -> %.3f deg with it",
              amp[0][0], amp[0][1], amp[1][0], amp[1][1]));
}

int main()
{
    test_torque_law();
    test_filter_pi();
    test_gains();
    test_tuning();
    test_closed_loop();
    test_input();
    test_floating();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
