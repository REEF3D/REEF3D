// REEF3D::CORAL unit test of the wind turbine rotor in the multibody model (coral/aero/coral_rotor, coral_rotor_aero,
// coral/structure/coral_drivetrain, coral_blade_table): a three-bladed rotor of the Glauert optimum design in
// uniform axial wind = the steady blade element momentum (rotor_steady) at any azimuth; cone (Vx and Vy of the coned
// element); spin-up of the free rotor against J domega/dt = Q(omega) of the steady rotor; the k omega^2 generator
// settles where Q(omega) = K omega^2; constant torque and the slip generator; torsion of the flexible drivetrain at
// its natural frequency; yawed inflow with the skewed wake (yaw angle, larger induction downwind); tower potential
// flow; mass and inertia of hub and blades; ElastoDyn's blade file; the keyword with its channels.
// Architect: Hans Bihs
#include"coral/aero/coral_rotor.h"
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

static rotor_def steady_def(const rotor_spec &sp)
{
    rotor_def rd;
    rd.B = sp.B; rd.R_hub = sp.R_hub; rd.R_tip = sp.stations.back().r;
    rd.stations = sp.stations; rd.airfoils = sp.airfoils;
    return rd;
}

static void test_steady()
{
    std::cout<<"rotor in uniform axial wind = steady BEM"<<std::endl;
    rotor_spec sp = design();
    auto w = std::make_shared<wind_field>(U0, 90.0);
    model m;
    m.gravity = Vec3::Zero();
    const rotor_handles h = build_rotor(m, w, sp);
    const rotor_result rs = rotor_steady(steady_def(sp), U0, OM, 0.0, 1.225);
    double e = 0.0;
    for(double az : {0.0, 17.0, 60.0, 133.0})
    {
        // turn the rotor body to the azimuth
        const Vec3 ex(1, 0, 0);
        m.bodies[h.hub].R = Eigen::AngleAxisd(az*deg, ex).toRotationMatrix();
        m.bodies[h.hub].w = OM*m.bodies[h.hub].R.transpose()*ex;
        std::vector<station_load> L;
        h.aero->station_loads(m, 0.0, L);
        const auto v = h.aero->values(m, 0.0, L);
        e = std::max(e, std::max(std::fabs(v.T/rs.T - 1.0), std::fabs(v.Q/rs.Q - 1.0)));
    }
    check(e < 1e-10, fmt("thrust %.1f kN, torque %.1f kN m (Cp %.4f) = rotor_steady at four azimuths (err %.1e)", rs.T/1e3, rs.Q/1e3, rs.Cp, e));
    check(std::fabs(h.R_tip - 40.0) < 1e-12 && h.aero->stations().size()==93, "93 stations on three blades, R_tip 40 m");
    // mass and inertia: hub + three blades of 38 m with 400 - 250 - 100 kg/m (trapezoid on the 3 stations)
    const double mb = 0.5*19.0*(400 + 250) + 0.5*19.0*(250 + 100);
    const double Jb = 0.5*19.0*400*4.0 + 19.0*250*21.0*21.0 + 0.5*19.0*100*40.0*40.0;
    check(std::fabs(h.mass - (50000.0 + 3*mb)) < 1e-6 && std::fabs(h.J - (1e5 + 3*Jb)) < 1e-3,
          fmt("rotor mass %.0f kg, inertia %.4e kg m^2 (hub + blade point masses)", h.mass, h.J));
}

static void test_cone()
{
    std::cout<<"cone and tilt"<<std::endl;
    rotor_spec sp = design();
    sp.precone = -2.5*deg;
    const double tilt = 5.0*deg;
    sp.axis = Vec3(std::cos(tilt), 0, -std::sin(tilt));
    auto w = std::make_shared<wind_field>(U0, 90.0);
    model m;
    const rotor_handles h = build_rotor(m, w, sp);
    std::vector<station_load> L;
    h.aero->station_loads(m, 0.0, L);
    // a coned element in tilted wind: Vx = U (cone/tilt), Vy = Omega r cos(cone) +- the in-plane wind, r in the plane
    const aero_station &s = h.aero->stations()[20];
    const Vec3 P = m.body(h.hub).point(s.s), n = m.body(h.hub).R*s.n0, t = m.body(h.hub).R*s.t0;
    const Vec3 V = Vec3(U0, 0, 0) - m.body(h.hub).point_velocity(s.s);
    const double rr = sp.stations[20].r*std::cos(sp.precone);
    const Vec3 d = P - sp.hub, e = sp.axis.normalized();
    check(std::fabs((d - d.dot(e)*e).norm() - rr) < 1e-9 && std::fabs(n.dot(e) - std::cos(sp.precone)) < 1e-12 && std::fabs(t.dot(e)) < 1e-12,
          fmt("blade 1 coned upwind by 2.5 deg: in-plane radius %.4f m, normal tilted by the cone", rr));
    bem_element el;
    el.r = rr; el.chord = s.chord; el.theta = s.twist; el.R = 40*std::cos(sp.precone); el.R_hub = 2*std::cos(sp.precone); el.B = 3;
    el.af = &sp.airfoils[0];
    bem_options o;
    // the skewed wake from the tilt (5 deg): its factor
    const bem_result b0 = bem_solve(el, V.dot(n), -V.dot(t), 1.225, o);
    const double gam = tilt, chi = (0.6*b0.a + 1)*gam;
    const Vec3 dr = d - d.dot(e)*e, cr = (Vec3(U0, 0, 0) - U0*std::cos(tilt)*e).normalized();
    const double K = 15*M_PI/32*(rr/el.R)*std::tan(0.5*chi)*dr.dot(cr)/rr;
    const bem_result b = bem_solve(el, V.dot(n), -V.dot(t), 1.225, o, K);
    check(std::fabs(L[20].bem.fn - b.fn) < 1e-9*std::fabs(b.fn) && std::fabs(L[20].bem.a - b.a) < 1e-12,
          fmt("the station's element = BEM of its own Vx %.3f, Vy %.3f m/s with the tilt's skew factor %.4f", V.dot(n), -V.dot(t), K));
    const auto v = h.aero->values(m, 0.0, L);
    check(std::fabs(v.yaw - tilt) < 1e-12, fmt("skew angle of the inflow = the tilt (%.2f deg)", v.yaw/deg));
}

static void test_spinup()
{
    std::cout<<"free rotor spin-up, k omega^2 generator, constant torque, slip"<<std::endl;
    rotor_spec sp = design(15);
    sp.omega0 = 0.6;
    auto w = std::make_shared<wind_field>(U0, 90.0);
    model m;
    m.gravity = Vec3::Zero();
    rotor_handles h = build_rotor(m, w, sp);
    timeint ti;
    ti.initialize(m, 0.0);
    const double dt = 0.05;
    // reference: J domega/dt = Q(omega) with RK4 on rotor_steady
    const rotor_def rd = steady_def(sp);
    auto Q = [&](double om) {return rotor_steady(rd, U0, om, 0.0, 1.225).Q;};
    double om = 0.6, emax = 0.0;
    for(int i=0; i<200; ++i)
    {
        ti.step(m, i*dt, dt);
        const double k1 = Q(om)/h.J, k2 = Q(om + 0.5*dt*k1)/h.J, k3 = Q(om + 0.5*dt*k2)/h.J, k4 = Q(om + dt*k3)/h.J;
        om += dt/6*(k1 + 2*k2 + 2*k3 + k4);
        emax = std::max(emax, std::fabs(h.drive->omega_rotor(m) - om));
    }
    check(emax < 1e-3*om, fmt("free rotor from 0.6 rad/s: %.4f rad/s after 10 s, the torque balance %.4f (largest difference %.1e)", h.drive->omega_rotor(m), om, emax));

    // k omega^2: equilibrium Q(omega) = K omega^2 with K = Q(OM)/OM^2
    rotor_spec s2 = design(15);
    const double K = Q(OM)/(OM*OM);
    s2.omega0 = 1.6;
    s2.dt.law = generator_law::k_omega2; s2.dt.K = K;
    model m2;
    m2.gravity = Vec3::Zero();
    rotor_handles h2 = build_rotor(m2, w, s2);
    timeint t2;
    t2.initialize(m2, 0.0);
    for(int i=0; i<3000; ++i) t2.step(m2, i*dt, dt);
    const double ws = h2.drive->omega_rotor(m2);
    check(std::fabs(ws/OM - 1.0) < 1e-4, fmt("k omega^2 generator settles at %.5f rad/s, where Q = K omega^2 (%.5f)", ws, OM));
    std::vector<std::string> nm; std::vector<double> vv;
    h2.drive->channel_names(nm); h2.drive->channel_values(m2, 150.0, vv);
    check(nm[3]=="rotor_dt_P_gen" && std::fabs(vv[3]/(Q(ws)*ws) - 1.0) < 1e-3, fmt("generator power %.1f kW = aerodynamic power", vv[3]/1e3));

    // gearbox 50: generator torque Q/50 on the high-speed shaft at 50 OM holds the speed (constant torque)
    rotor_spec s3 = design(15);
    s3.dt.N = 50.0; s3.dt.law = generator_law::constant; s3.dt.T = Q(OM)/50.0; s3.J_gen = 20.0;
    model m3;
    m3.gravity = Vec3::Zero();
    rotor_handles h3 = build_rotor(m3, w, s3);
    timeint t3;
    t3.initialize(m3, 0.0);
    for(int i=0; i<100; ++i) t3.step(m3, i*dt, dt);
    check(std::fabs(h3.drive->omega_rotor(m3)/OM - 1.0) < 1e-6 && std::fabs(h3.drive->omega_hss(m3)/(50*OM) - 1.0) < 1e-6,
          "constant torque Q/N on the high-speed shaft (N = 50) holds the design speed");
    // slip generator with a steep slope: the speed stays near the synchronous one
    rotor_spec s4 = design(15);
    s4.omega0 = 1.0;
    s4.dt.law = generator_law::slip; s4.dt.slope = 1.0e8; s4.dt.omega_sync = 1.2;
    model m4;
    m4.gravity = Vec3::Zero();
    rotor_handles h4 = build_rotor(m4, w, s4);
    timeint t4;
    t4.initialize(m4, 0.0);
    for(int i=0; i<200; ++i) t4.step(m4, i*dt, dt);
    const double sl = h4.drive->omega_rotor(m4) - 1.2;
    check(std::fabs(sl - Q(1.2 + sl)/1.0e8) < 1e-9, fmt("slip generator: slip %.2e rad/s = Q/slope", sl));
}

static void test_torsion()
{
    std::cout<<"flexible drivetrain"<<std::endl;
    rotor_spec sp = design(15);
    sp.omega0 = 0.0;
    sp.aero.rho = 0.0;                 // no aerodynamic load
    sp.dt.N = 97.0; sp.J_gen = 534.116; sp.dt.k = 8.676e8; sp.dt.c = 0.0;
    auto w = std::make_shared<wind_field>(0.0, 90.0);
    model m;
    m.gravity = Vec3::Zero();
    rotor_handles h = build_rotor(m, w, sp);
    m.bodies[h.gen].w = Vec3(0.01, 0, 0);     // the generator side starts to turn
    timeint_settings ts;
    ts.rho_inf = 1.0;
    timeint ti(ts);
    ti.initialize(m, 0.0);
    const double Jr = h.J, Jg = 534.116*97*97;
    const double fn = std::sqrt(8.676e8*(Jr + Jg)/(Jr*Jg))/(2*M_PI);
    const double dt = 1.0/(fn*200.0);
    double tp = 0.0, last = 0.0, prev = h.drive->shaft_torque(m);
    int nz = 0;
    for(int i=0; i<2000; ++i)
    {
        ti.step(m, i*dt, dt);
        const double T = h.drive->shaft_torque(m);
        if(prev < 0.0 && T >= 0.0)
        {
            const double tz = (i + prev/(prev - T))*dt;
            if(nz++ > 0) tp = tz - last;
            last = tz;
        }
        prev = T;
    }
    check(nz > 3 && std::fabs(1.0/tp/fn - 1.0) < 2e-3, fmt("torsion at %.4f Hz = sqrt(k (J_r + N^2 J_g)/(J_r N^2 J_g))/2 pi = %.4f Hz", 1.0/tp, fn));
    const double Lr = Jr*h.drive->omega_rotor(m) + Jg*m.body(h.gen).w(0);
    check(std::fabs(Lr - Jg*0.01) < 1e-6*Jg*0.01, "angular momentum of rotor and generator kept");
}

static void test_yaw_tower()
{
    std::cout<<"yawed inflow, tower potential flow"<<std::endl;
    rotor_spec sp = design();
    auto w = std::make_shared<wind_field>(U0, 90.0, 30.0*deg);
    model m;
    const rotor_handles h = build_rotor(m, w, sp);
    std::vector<station_load> L;
    // blade 1 horizontal: on the downwind side of the cross flow (the cross flow points to +y, blade 1 at azimuth
    // -90 deg points along +y for the clockwise rotor)
    m.bodies[h.hub].R = Eigen::AngleAxisd(-90.0*deg, Vec3::UnitX()).toRotationMatrix();
    h.aero->station_loads(m, 0.0, L);
    const auto v = h.aero->values(m, 0.0, L);
    const Vec3 P1 = L[20].P - sp.hub;
    // the same station on the opposite side: rotate the rotor by 180 deg
    m.bodies[h.hub].R = Eigen::AngleAxisd(90.0*deg, Vec3::UnitX()).toRotationMatrix();
    std::vector<station_load> L2;
    h.aero->station_loads(m, 0.0, L2);
    check(std::fabs(v.yaw - 30*deg) < 1e-12 && P1(1) > 0.0 && L[20].bem.a > L2[20].bem.a,
          fmt("yaw 30 deg: induction a %.4f on the downwind side, %.4f upwind", L[20].bem.a, L2[20].bem.a));
    const rotor_result rs = rotor_steady(steady_def(sp), U0, OM, 0.0, 1.225);
    check(v.T < rs.T && v.T > 0.6*rs.T, fmt("thrust %.0f kN, %.0f kN without yaw", v.T/1e3, rs.T/1e3));

    rotor_spec st = design();
    st.aero.tower = true;
    st.aero.tower_base = Vec3(5.0, 0, 0); st.aero.tower_top = Vec3(5.0, 0, 88.0);
    st.aero.D_base = 6.0; st.aero.D_top = 4.0;
    auto w0 = std::make_shared<wind_field>(U0, 90.0);
    model mt;
    const rotor_handles ht = build_rotor(mt, w0, st);
    const Vec3 u(U0, 0, 0);
    // at z = 44 the diameter is 5 m: 2 radii upstream and 2 radii to the side
    const Vec3 a = ht.aero->tower_influence(mt, Vec3(0.0, 0, 44.0), u), b = ht.aero->tower_influence(mt, Vec3(5.0, 5.0, 44.0), u);
    const Vec3 c = ht.aero->tower_influence(mt, Vec3(2.5, 2.5, 44.0), u), d = ht.aero->tower_influence(mt, Vec3(0.0, 0, 95.0), u);
    check(std::fabs(a(0) + 0.25*U0) < 1e-12 && std::fabs(b(0) - 0.25*U0) < 1e-12 && std::fabs(b(1)) < 1e-12
          && std::fabs(c(1) + 2*U0*(-1.0*1.0)/4.0) < 1e-12 && d.norm()==0.0,
          "potential flow: u 0.75 U two radii upstream, 1.25 U beside, v at 45 deg, nothing above the tower top");
}

static void test_input()
{
    std::cout<<"keyword"<<std::endl;
    {
        std::ofstream f("rotor_test_af.dat");
        f<<"! plain polar\n";
        const airfoil af = thin();
        for(size_t i=0; i<af.tables[0].alpha.size(); ++i)
            f<<af.tables[0].alpha[i]/deg<<" "<<af.tables[0].cl[i]<<" "<<af.tables[0].cd[i]<<" "<<af.tables[0].cm[i]<<"\n";
        std::ofstream b("rotor_test_blade.dat");
        const rotor_spec sp = design(10);
        b<<"------- AERODYN v15.00.* BLADE DEFINITION INPUT FILE ----\ntest\n====== Blade Properties ====\n"
         <<sp.stations.size()<<"   NumBlNds\nBlSpn BlCrvAC BlSwpAC BlCrvAng BlTwist BlChord BlAFID\n(m) (m) (m) (deg) (deg) (m) (-)\n";
        for(const blade_station &s : sp.stations)
            b<<s.r - 2.0<<" 0 0 0 "<<s.twist/deg<<" "<<s.chord<<" 1\n";
        std::ofstream e("rotor_test_ed.dat");
        e<<"------- ELASTODYN V1.00.* INDIVIDUAL BLADE INPUT FILE ----\ntest\n---- BLADE PARAMETERS ----\n"
         <<"          3   NBlInpSt    - Number of blade input stations (-)\n     1.0   AdjBlMs\n     1.0   AdjFlSt\n     1.0   AdjEdSt\n"
         <<"---- DISTRIBUTED BLADE PROPERTIES ----\nBlFract PitchAxis StrcTwst BMassDen FlpStff EdgStff\n(-) (-) (deg) (kg/m) (Nm^2) (Nm^2)\n"
         <<"0.0 0.25 13.3 400 1e10 1e10\n0.5 0.30 5.0 250 1e9 2e9\n1.0 0.25 0.1 100 1e7 2e7\n"
         <<"---- BLADE MODE SHAPES ----\n    0.0622   BldFl1Sh(2)\n    1.7254   BldFl1Sh(3)\n";
    }
    coral_case c;
    c.read_string("time dt 0.05 end 1\n"
                  "gravity 0 0 0\n"
                  "wind steady U 10 z_ref 90\n"
                  "body nacelle mass 1000 cog 0 0 90\n"
                  "joint yaw fixed ground nacelle point 0 0 90\n"
                  "rotor rotor hub 0 0 90 axis 1 0 0 R_hub 2 blade rotor_test_blade.dat airfoils 1 rotor_test_af.dat \\\n"
                  "   structure rotor_test_ed.dat hub_mass 50000 hub_inertia 1e5 nacelle nacelle rpm 16.71127 \\\n"
                  "   gearbox 50 gen_inertia 20 generator k_omega2 0.1\n", ".");
    std::remove("rotor_test_af.dat"); std::remove("rotor_test_blade.dat"); std::remove("rotor_test_ed.dat");
    load *ra = nullptr, *dt = nullptr;
    for(auto &l : c.m.load_list)
    {
        if(l->name()=="rotor") ra = l.get();
        if(l->name()=="rotor_dt") dt = l.get();
    }
    check(ra && dt && c.m.find_body("rotor_hub") >= 0 && c.joints.count("rotor_shaft"), "rotor, drivetrain, hub body and shaft joint from the keyword");
    if(!ra) return;
    std::vector<std::string> n;
    std::vector<double> v;
    ra->channel_names(n);
    ra->channel_values(c.m, 0.0, v);
    const rotor_spec sp = design(10);
    const rotor_result rs = rotor_steady(steady_def(sp), U0, OM, 0.0, 1.225);
    check(n.size()==17 && n[2]=="rotor_T" && std::fabs(v[1] - OM*30/M_PI) < 1e-3 && std::fabs(v[2]/rs.T - 1.0) < 1e-5 && std::fabs(v[5]/rs.Cp - 1.0) < 1e-5,
          fmt("channels: rpm %.3f, thrust %.1f kN, Cp %.4f as the steady rotor", v[1], v[2]/1e3, v[5]));
    check(v[11] > 0.0 && v[12] > 0.0, fmt("blade 1 root moments out of plane %.0f kN m, in plane %.0f kN m", v[11]/1e3, v[12]/1e3));
    auto err = [](const std::string &s)
    {
        try {coral_case k; k.read_string(s, ".");} catch(const std::exception &e) {return std::string(e.what());}
        return std::string();
    };
    check(err("rotor r hub 0 0 90 axis 1 0 0 R_hub 2 blade x.dat airfoils 1 a.dat\n").find("needs a wind") != std::string::npos, "rotor without wind: error");
}

int main()
{
    test_steady();
    test_cone();
    test_spinup();
    test_torsion();
    test_yaw_tower();
    test_input();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
