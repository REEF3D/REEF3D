// REEF3D::CORAL unit test: unsteady aerodynamics (src/coral/aero/coral_unsteady_aero, coral_rotor_aero): the exact
// lag step, Øye's dynamic inflow (step and frequency response against the closed forms), the static polar prepared
// for dynamic stall, Øye's and the Beddoes-Leishman (HGM) dynamic stall model (the steady state = the static polar,
// the attached-flow response = the indicial function of the model and Theodorsen's function, Øye's step response, the stall
// hysteresis and the quasi-steady limit), and the rotor with dynamic inflow and dynamic stall in the multibody model
// (a steady rotor stays steady, the induced velocity after a pitch step = the closed form of the filter, no load at
// the tip station wherever it moves, restart).
// Architect: Hans Bihs

#include"coral/aero/coral_unsteady_aero.h"
#include"coral/aero/coral_rotor.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/io/coral_input.h"
#include"coral/wind/coral_wind.h"
#include<cmath>
#include<complex>
#include<cstdio>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>

using namespace coral;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  ok    " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}
static std::string fmt(const char *f, double a, double b=0.0, double c=0.0, double d=0.0)
{
    char s[512];
    std::snprintf(s, sizeof(s), f, a, b, c, d);
    return s;
}

static const double deg = M_PI/180.0;

// 2 pi alpha up to 12 deg, a flat plate beyond (as coral_rotor_test)
static airfoil thin()
{
    airfoil af;
    af.tables.emplace_back();
    polar_table &p = af.tables.back();
    for(int i=-180; i<=180; ++i)
    {
        const double a = i*deg;
        const double cl = std::fabs(a) < 12*deg ? 2*M_PI*a : 2*M_PI*12*deg*(a > 0 ? 1 : -1)*std::cos(a)*std::cos(a)/std::pow(std::cos(12*deg), 2);
        p.alpha.push_back(a); p.cl.push_back(cl); p.cd.push_back(0.008 + 1.2*std::sin(a)*std::sin(a)); p.cm.push_back(-0.05*std::sin(a));
    }
    return af;
}

static void test_lag()
{
    std::cout<<"exact lag step"<<std::endl;
    // y' = (u - y)/tau with u = 1 + 2 t, y(0) = 3: y = u - 2 tau + (3 - 1 + 2 tau) exp(-t/tau)
    const double tau = 0.7;
    double e = 0.0;
    for(double dt : {0.001, 0.1, 0.5, 3.0})
    {
        double y = 3.0, t = 0.0;
        for(int i=0; i<int(6.0/dt); ++i)
        {
            y = lag_step(y, 1 + 2*t, 1 + 2*(t + dt), tau, dt);
            t += dt;
            e = std::max(e, std::fabs(y - (1 + 2*t - 2*tau + (2 + 2*tau)*std::exp(-t/tau))));
        }
    }
    check(e < 1e-12, fmt("linear input, dt/tau from 0.0014 to 4.3: exact (largest error %.1e)", e));
    check(lag_step(5.0, 1.0, 2.0, 1.0, 0.0)==5.0 && std::fabs(lag_step(5.0, 1.0, 2.0, 1.0e-9, 1.0) - 2.0) < 1e-8,
          "dt = 0 keeps the state, tau -> 0 gives the input");
}

static void test_inflow()
{
    std::cout<<"Øye's dynamic inflow"<<std::endl;
    check(std::fabs(oye_tau1(0.3, 60.0, 10.0) - 1.1/(1 - 0.39)*6.0) < 1e-12 && std::fabs(oye_tau1(0.7, 60.0, 10.0) - 1.1/0.35*6.0) < 1e-12
          && std::fabs(oye_tau2(2.0, 0.5) - (0.39 - 0.065)*2.0) < 1e-12, "tau1 = 1.1/(1 - 1.3 a) R/U (a up to 0.5), tau2 = (0.39 - 0.26 (r/R)^2) tau1");
    // step of W_qs from 0 to 1: W_int = 1 - (1 - k) exp(-t/tau1),
    // W = 1 - exp(-t/tau2) - (1 - k) tau1/(tau1 - tau2) (exp(-t/tau1) - exp(-t/tau2))
    const double t1 = 4.0, t2 = 1.2, k = 0.6, dt = 1.0e-4;
    const double z[2] = {0.0, 0.0}, one[2] = {1.0, 0.5};
    inflow_state s = inflow_steady(z);
    double e = 0.0, t = 0.0;
    for(int i=0; i<int(20.0/dt); ++i)
    {
        s = inflow_step(s, one, t1, t2, k, dt);
        t += dt;
        const double Wi = 1 - (1 - k)*std::exp(-t/t1);
        const double W = 1 - std::exp(-t/t2) - (1 - k)*t1/(t1 - t2)*(std::exp(-t/t1) - std::exp(-t/t2));
        e = std::max({e, std::fabs(s.Wint[0] - Wi), std::fabs(s.W[0] - W), std::fabs(s.W[1] - 0.5*W)});
    }
    check(e < 2e-4, fmt("step of W_qs: W_int and W = the closed forms (largest error %.1e, the step spread over dt = 1e-4 s)", e));
    // harmonic W_qs: W/W_qs = (1 + k tau1 i w)/((1 + tau1 i w)(1 + tau2 i w))
    for(double w : {0.1, 0.5, 2.0})
    {
        const double h = 0.01;
        inflow_state q = inflow_steady(z);
        std::complex<double> acc = 0.0;
        const int n0 = int(10*t1/h);
        const double T = 2*M_PI/w;
        int cnt = 0;
        double tt = 0.0;
        for(int i=0; tt < 10*t1 + 10*T; ++i)
        {
            const double in[2] = {std::sin(w*(tt + h)), 0.0};
            q = inflow_step(q, in, t1, t2, k, h);
            tt += h;
            if(i >= n0 && tt <= 10*t1 + 10*T)
            {
                acc += q.W[0]*std::exp(std::complex<double>(0, -w*tt));
                ++cnt;
            }
        }
        const std::complex<double> got = acc*2.0/double(cnt)*std::complex<double>(0, 1);
        const std::complex<double> iw(0, w);
        const std::complex<double> ex = (1.0 + k*t1*iw)/((1.0 + t1*iw)*(1.0 + t2*iw));
        check(std::abs(got - ex) < 3e-3, fmt("harmonic W_qs at %.1f rad/s: W/W_qs = %.4f at %.1f deg (the transfer function: %.4f",
              w, std::abs(got), std::arg(got)/deg, std::abs(ex)) + fmt(" at %.1f deg)", std::arg(ex)/deg));
    }
}

static void test_polar()
{
    std::cout<<"the static polar for dynamic stall"<<std::endl;
    const airfoil af = thin();
    dynamic_stall_settings s;
    s.model = 2;
    const stall_polar p(af, 0.0, s);
    check(std::fabs(p.cl_alpha - 2*M_PI) < 1e-9 && std::fabs(p.alpha0) < 1e-12 && p.A1==0.165 && p.A2==0.335 && p.b1==0.0455 && p.b2==0.3
          && p.T_f0==3.0 && p.T_p==1.7, fmt("Cl_alpha = 2 pi (%.6f), alpha0 = 0, Jones' constants, T_f0 3, T_p 1.7", p.cl_alpha));
    double e = 0.0;
    bool f1 = true, mono = true;
    double fprev = 1.0;
    for(int i=-1800; i<=1800; ++i)
    {
        const double a = i*0.1*deg + 0.0123*deg;
        const double f = p.f_st(a);
        e = std::max(e, std::fabs(f*p.cl_inv(a) + (1 - f)*p.cl_fs(a) - af.coefficients(a).cl));
        if(std::fabs(a) <= 12*deg && std::fabs(f - 1.0) > 1e-12) f1 = false;
        if(a > 12*deg) {if(f > fprev + 1e-15) mono = false; fprev = f;}
    }
    check(e < 1e-12, fmt("f_st Cl_inv + (1 - f_st) Cl_fs = the static cl at every angle, between the table's angles too (%.1e)", e));
    check(f1 && mono && p.f_st(90*deg)==0.0 && p.f_st(20*deg) > 0.0 && p.f_st(20*deg) < 1.0,
          fmt("f_st = 1 in the linear range, falls beyond 12 deg (%.3f at 20 deg), 0 at 90 deg", p.f_st(20*deg)));
    // the airfoil file's constants
    airfoil b = af;
    b.tables[0].params["A1"] = 0.2; b.tables[0].params["T_f0"] = 4.0; b.tables[0].params["alpha0"] = -1.0;
    dynamic_stall_settings s2 = s;
    s2.T_f0 = 5.0;
    const stall_polar q(b, 0.0, s2);
    check(q.A1==0.2 && q.T_f0==5.0 && std::fabs(q.alpha0 + deg) < 1e-12, "A1 and alpha0 from the airfoil file, T_f0 of the input over the file's");
    dynamic_stall_settings s3;
    s3.model = 1;
    check(stall_polar(af, 0.0, s3).T_f0==6.0, "Øye's model: T_f0 6");
}

static void test_stall()
{
    std::cout<<"dynamic stall"<<std::endl;
    const airfoil af = thin();
    for(int model : {1, 2})
    {
        dynamic_stall_settings s;
        s.model = model;
        const stall_polar p(af, 0.0, s);
        double e = 0.0;
        for(int i=-60; i<=60; ++i)
        {
            stall_input in;
            in.alpha = i*0.5*deg + 0.031*deg; in.U = 50.0; in.chord = 2.0;
            stall_state x = stall_steady(p, model, in), y;
            airfoil_coef c{};
            for(int j=0; j<5; ++j) {c = stall_step(p, model, x, in, 0.01, y); x = y;}
            const airfoil_coef st = af.coefficients(in.alpha);
            e = std::max({e, std::fabs(c.cl - st.cl), std::fabs(c.cd - st.cd), std::fabs(c.cm - st.cm)});
        }
        check(e < 1e-12, fmt(model==1 ? "Øye: the steady state gives the static polar from -30 to 30 deg (%.1e)"
                                       : "HGM: the steady state gives the static cl, cd, cm from -30 to 30 deg (%.1e)", e));
    }
    // HGM, attached flow (2 +- 1 deg, no pitch rate): cl/(Cl_alpha alpha) = 1 - A1 - A2 + A1 b1/(b1 + i k) + A2 b2/(b2 + i k)
    {
        dynamic_stall_settings s;
        s.model = 2;
        const stall_polar p(af, 0.0, s);
        const double U = 40.0, c = 2.0, Tu = c/(2*U);
        for(double kr : {0.05, 0.2, 0.5})
        {
            const double w = kr/Tu, T = 2*M_PI/w, dt = Tu/50;
            stall_input in;
            in.U = U; in.chord = c; in.alpha = 2*deg;
            stall_state x = stall_steady(p, 2, in), y;
            std::complex<double> acc = 0.0;
            int cnt = 0;
            const int n = int(std::round(40*T/dt));
            for(int i=1; i<=n; ++i)
            {
                const double t = i*dt;
                in.alpha = 2*deg + 1*deg*std::sin(w*t);
                const airfoil_coef cc = stall_step(p, 2, x, in, dt, y);
                x = y;
                if(t > 20*T) {acc += (cc.cl - 2*M_PI*2*deg)*std::exp(std::complex<double>(0, -w*t)); ++cnt;}
            }
            const std::complex<double> got = acc*2.0/double(cnt)*std::complex<double>(0, 1)/(2*M_PI*1*deg);
            const std::complex<double> ik(0, kr);
            const std::complex<double> ex = 1.0 - p.A1 - p.A2 + p.A1*p.b1/(p.b1 + ik) + p.A2*p.b2/(p.b2 + ik);
            // Theodorsen's function C(k) = H1(k)/(H1(k) + i H0(k)), H = J - i Y (Hankel functions of the second kind)
            const std::complex<double> H1(std::cyl_bessel_j(1.0, kr), -std::cyl_neumann(1.0, kr));
            const std::complex<double> H0(std::cyl_bessel_j(0.0, kr), -std::cyl_neumann(0.0, kr));
            const std::complex<double> Ck = H1/(H1 + std::complex<double>(0, 1)*H0);
            check(std::abs(got - ex) < 2e-3 && std::abs(got - Ck) < 0.02,
                  fmt("HGM attached flow at k = %.2f: cl/(Cl_alpha alpha) = %.4f at %.2f deg (the model's indicial function: %.4f",
                  kr, std::abs(got), std::arg(got)/deg, std::abs(ex)) + fmt(" at %.2f deg; Theodorsen: %.4f at %.2f deg)", std::arg(ex)/deg, std::abs(Ck), std::arg(Ck)/deg));
        }
    }
    // Øye: step 5 -> 20 deg: f = f_st(20) + (1 - f_st(20)) exp(-t/T_f), cl = f Cl_inv + (1 - f) Cl_fs
    {
        dynamic_stall_settings s;
        s.model = 1;
        const stall_polar p(af, 0.0, s);
        const double U = 40.0, c = 2.0, Tf = p.T_f0*c/(2*U), dt = Tf*1e-3;
        stall_input in;
        in.U = U; in.chord = c; in.alpha = 5*deg;
        stall_state x = stall_steady(p, 1, in), y;
        in.alpha = 20*deg;
        const double fs = p.f_st(20*deg);
        double e = 0.0, cl0 = 0.0;
        for(int i=1; i<=int(5*Tf/dt); ++i)
        {
            const airfoil_coef cc = stall_step(p, 1, x, in, dt, y);
            x = y;
            if(i==1) cl0 = cc.cl;
            const double f = fs + (1 - fs)*std::exp(-i*dt/Tf);
            e = std::max(e, std::fabs(cc.cl - (f*p.cl_inv(20*deg) + (1 - f)*p.cl_fs(20*deg))));
        }
        check(e < 2e-3, fmt("Øye: step 5 -> 20 deg: cl from %.3f (attached) to the static %.3f in T_f, the closed form to %.1e",
              cl0, af.coefficients(20*deg).cl, e));
    }
    // stall hysteresis (15 +- 10 deg, k = 0.1) and the quasi-steady limit (k = 0.001)
    for(int model : {1, 2})
    {
        dynamic_stall_settings s;
        s.model = model;
        const stall_polar p(af, 0.0, s);
        const double U = 40.0, c = 2.0, Tu = c/(2*U);
        double clmax[2] = {0, 0}, dqs = 0.0, up = 0.0, down = 0.0;
        for(int q=0; q<2; ++q)
        {
            const double kr = q==0 ? 0.1 : 0.001, w = kr/Tu, T = 2*M_PI/w, dt = std::min(Tu/20, T/2000);
            stall_input in;
            in.U = U; in.chord = c; in.alpha = 15*deg;
            stall_state x = stall_steady(p, model, in), y;
            const int n = int(std::round(3*T/dt));
            for(int i=1; i<=n; ++i)
            {
                const double t = i*dt;
                in.alpha = 15*deg + 10*deg*std::sin(w*t);
                const airfoil_coef cc = stall_step(p, model, x, in, dt, y);
                x = y;
                if(t < 2*T) continue;
                clmax[q] = std::max(clmax[q], cc.cl);
                if(q==1) dqs = std::max(dqs, std::fabs(cc.cl - af.coefficients(in.alpha).cl));
                if(q==0 && std::fabs(std::fmod(t, T)) < dt) up = cc.cl;                          // 15 deg, upstroke
                if(q==0 && std::fabs(std::fmod(t, T) - 0.5*T) < 0.5*dt) down = cc.cl;           // 15 deg, downstroke
            }
        }
        const double stmax = af.coefficients(12*deg).cl;
        check(clmax[0] > 1.1*stmax && up > down + 0.1 && dqs < 0.03,
              fmt(model==1 ? "Øye, 15 +- 10 deg at k = 0.1: cl max %.3f (static %.3f), at 15 deg up %.3f, down %.3f;"
                           : "HGM, 15 +- 10 deg at k = 0.1: cl max %.3f (static %.3f), at 15 deg up %.3f, down %.3f;", clmax[0], stmax, up, down)
              + fmt(" k = 0.001: the static polar to %.3f", dqs));
    }
}

// the Glauert-optimum rotor of coral_rotor_test (R 40 m, tip speed ratio 7 at 10 m/s), at a constant speed
static rotor_spec design(int n=20)
{
    rotor_spec sp;
    sp.B = 3;
    sp.hub = Vec3(0, 0, 90.0);
    sp.axis = Vec3(1, 0, 0);
    sp.R_hub = 2.0;
    sp.airfoils.push_back(thin());
    const double ad = 6.0*deg, cld = 2*M_PI*ad, R = 40.0, LAM = 7.0;
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
    sp.hub_J = 1.0e14;                  // the speed stays constant
    sp.omega0 = LAM*10.0/R;
    sp.aero.rho = 1.225;
    return sp;
}

static void test_rotor()
{
    std::cout<<"rotor with dynamic inflow and dynamic stall"<<std::endl;
    auto w = std::make_shared<wind_field>(10.0, 90.0);
    // steady: with both models the loads stay the quasi-steady ones
    {
        rotor_spec a = design(), b = design();
        b.aero.dyn_inflow.on = true;
        b.aero.dyn_stall.model = 2;
        model ma, mb;
        ma.gravity = mb.gravity = Vec3::Zero();
        const rotor_handles ha = build_rotor(ma, w, a), hb = build_rotor(mb, w, b);
        timeint ta, tb;
        ta.initialize(ma, 0.0); tb.initialize(mb, 0.0);
        const double dt = 0.02;
        double e = 0.0;
        for(int i=0; i<100; ++i)
        {
            ta.step(ma, i*dt, dt); tb.step(mb, i*dt, dt);
            std::vector<station_load> La, Lb;
            ha.aero->station_loads(ma, (i + 1)*dt, La); hb.aero->station_loads(mb, (i + 1)*dt, Lb);
            const auto va = ha.aero->values(ma, (i + 1)*dt, La), vb = hb.aero->values(mb, (i + 1)*dt, Lb);
            e = std::max({e, std::fabs(vb.T/va.T - 1), std::fabs(vb.Q/va.Q - 1)});
        }
        check(e < 1e-7 && hb.aero->tau1() > 0.0, fmt("steady inflow: thrust and torque = quasi-steady over 2 s (%.1e); tau1 = %.2f s", e, hb.aero->tau1()));
    }
    // pitch step 0 -> 4 deg over 0.1 s at 1 s: after it W_qs is constant and W follows the closed form of the filter
    {
        rotor_spec b = design();
        b.aero.dyn_inflow.on = true;
        b.pitch_t = {0.0, 1.0, 1.1, 100.0};
        b.pitch_v = {0.0, 0.0, 4*deg, 4*deg};
        rotor_spec q = b;
        q.aero.dyn_inflow.on = false;
        model mb, mq;
        mb.gravity = mq.gravity = Vec3::Zero();
        const rotor_handles hb = build_rotor(mb, w, b), hq = build_rotor(mq, w, q);
        timeint tb, tq;
        tb.initialize(mb, 0.0); tq.initialize(mq, 0.0);
        const double dt = 0.01;
        const int is = 12;
        inflow_state s1;
        double tau1 = 0.0, tau2 = 0.0, e = 0.0, D0 = 0.0, T0 = 0.0, Tq0 = 0.0, T11 = 0.0, Tq11 = 0.0, Tend = 0.0, Tqend = 0.0;
        double Wqs = 0.0;
        for(int i=0; i<4000; ++i)
        {
            const double t = (i + 1)*dt;
            tb.step(mb, i*dt, dt); tq.step(mq, i*dt, dt);
            std::vector<station_load> Lb, Lq;
            hb.aero->station_loads(mb, t, Lb); hq.aero->station_loads(mq, t, Lq);
            const double T = hb.aero->values(mb, t, Lb).T, Tq = hq.aero->values(mq, t, Lq).T;
            if(i==99) {T0 = T; Tq0 = Tq;}
            const inflow_state &s = hb.aero->inflow_states()[is];
            if(i==110)
            {
                s1 = s; tau1 = hb.aero->tau1();
                tau2 = oye_tau2(tau1, (b.stations[is].r)/40.0);
                Wqs = s.Wqs[0];
                D0 = s.Wint[0] - Wqs;
                T11 = T; Tq11 = Tq;
            }
            if(i > 110)
            {
                const double x = t - 1.11, E0 = s1.W[0] - Wqs;
                const double c1 = D0*tau1/(tau1 - tau2);
                const double W = Wqs + c1*std::exp(-x/tau1) + (E0 - c1)*std::exp(-x/tau2);
                e = std::max(e, std::fabs(s.W[0] - W)/std::fabs(s1.W[0] - Wqs));
            }
            Tend = T; Tqend = Tq;
        }
        check(e < 2e-3 && std::fabs(hb.aero->tau1()/tau1 - 1) < 1e-6,
              fmt("pitch step 0 -> 4 deg: the induced velocity relaxes as the closed form of the filter (tau1 %.2f s, tau2 %.2f s; %.1e of its change)", tau1, tau2, e));
        check(std::fabs(T0/Tq0 - 1) < 1e-7 && (Tq11 - Tq0) < 0.0 && (T11 - T0) < 1.2*(Tq11 - Tq0) && std::fabs(Tend/Tqend - 1) < 2e-3,
              fmt("thrust: %.1f kN before; after the step %.1f kN with dynamic inflow (frozen wake), %.1f kN quasi-steady;", T0/1e3, T11/1e3, Tq11/1e3)
              + fmt(" after 39 s %.1f = %.1f kN", Tend/1e3, Tqend/1e3));
    }
    // the tip station on a pitch bearing: no load wherever rounding or a small motion puts it (it switched between no
    // load and a load with a tiny tip loss factor, which stalled Newton)
    {
        rotor_spec b = design();
        b.pitch_joints = true;
        b.pitch = 8*deg;
        model mb;
        mb.gravity = Vec3::Zero();
        const rotor_handles hb = build_rotor(mb, w, b);
        const auto &st = hb.aero->stations();
        int tip = -1;
        for(size_t i=0; i<st.size(); ++i) if(st[i].blade==0) tip = (int)i;
        double fmax = 0.0;
        for(double d : {-1e-6, -1e-8, 0.0, 1e-8, 1e-6})
        {
            rigid_body B0 = mb.bodies[st[tip].body];
            mb.bodies[st[tip].body].x += d*(mb.bodies[st[tip].body].R*st[tip].z0);
            std::vector<station_load> L;
            hb.aero->station_loads(mb, 0.0, L);
            fmax = std::max(fmax, L[tip].F.norm());
            mb.bodies[st[tip].body] = B0;
        }
        check(fmax==0.0, "the tip station on a pitch bearing moved by -1e-6 .. 1e-6 m along the span: no load");
    }
    // restart: the state of the unsteady aerodynamics saved and loaded gives the same loads
    {
        rotor_spec b = design();
        b.aero.dyn_inflow.on = true;
        b.aero.dyn_stall.model = 2;
        b.pitch_t = {0.0, 0.5, 0.6, 100.0};
        b.pitch_v = {0.0, 0.0, 20*deg, 20*deg};
        model mb, mc;
        mb.gravity = mc.gravity = Vec3::Zero();
        const rotor_handles hb = build_rotor(mb, w, b), hc = build_rotor(mc, w, b);
        timeint tb;
        tb.initialize(mb, 0.0);
        for(int i=0; i<70; ++i) tb.step(mb, i*0.01, 0.01);
        std::vector<double> v;
        hb.aero->save_state(v);
        size_t k = 0;
        timeint tc;
        tc.initialize(mc, 0.0);
        const bool ok = hc.aero->load_state(v, k) && k==v.size();
        std::vector<station_load> L1, L2;
        hb.aero->station_loads(mb, 0.72, L1); hc.aero->station_loads(mb, 0.72, L2);
        double e = 0.0;
        for(size_t i=0; i<L1.size(); ++i) e = std::max(e, (L1[i].F - L2[i].F).norm());
        check(ok && e==0.0 && v.size()==4 + 14*L1.size(), fmt("restart: %.0f values saved, the same loads after loading", (double)v.size()));
    }
    // input
    {
        {
            std::ofstream f("ua_test_af.dat");
            const airfoil af = thin();
            for(size_t i=0; i<af.tables[0].alpha.size(); ++i)
                f<<af.tables[0].alpha[i]/deg<<" "<<af.tables[0].cl[i]<<" "<<af.tables[0].cd[i]<<" "<<af.tables[0].cm[i]<<"\n";
            std::ofstream b("ua_test_blade.dat");
            const rotor_spec sp = design(10);
            b<<"------- AERODYN v15.00.* BLADE DEFINITION INPUT FILE ----\ntest\n====== Blade Properties ====\n"
             <<sp.stations.size()<<"   NumBlNds\nBlSpn BlCrvAC BlSwpAC BlCrvAng BlTwist BlChord BlAFID\n(m) (m) (m) (deg) (deg) (m) (-)\n";
            for(const blade_station &x : sp.stations)
                b<<x.r - 2.0<<" 0 0 0 "<<x.twist/deg<<" "<<x.chord<<" 1\n";
        }
        const std::string head = "time dt 0.05 end 1\ngravity 0 0 0\nwind steady U 10 z_ref 90\nbody nacelle mass 1000 cog 0 0 90\n"
                                 "joint yaw fixed ground nacelle point 0 0 90\n"
                                 "rotor rotor hub 0 0 90 axis 1 0 0 R_hub 2 blade ua_test_blade.dat airfoils 1 ua_test_af.dat hub_mass 50000 "
                                 "hub_inertia 1e5 nacelle nacelle rpm 16.7";
        auto err = [](const std::string &s)
        {
            try {coral_case k; k.read_string(s, ".");} catch(const std::exception &e) {return std::string(e.what());}
            return std::string();
        };
        coral_case c;
        c.read_string(head + " dynamic_inflow di_k 0.5 dynamic_stall hgm ds_T_f0 4\n", ".");
        const rotor_aero *ra = c.rotors.count("rotor") ? c.rotors["rotor"].aero : nullptr;
        std::vector<std::string> n;
        if(ra) ra->channel_names(n);
        check(ra && ra->settings().dyn_inflow.on && ra->settings().dyn_inflow.k==0.5 && ra->settings().dyn_stall.model==2
              && ra->settings().dyn_stall.T_f0==4.0 && n.back()=="rotor_tau1",
              "keyword rotor ... dynamic_inflow di_k 0.5 dynamic_stall hgm ds_T_f0 4: the settings, channel rotor_tau1");
        coral_case d;
        d.read_string(head + "\n", ".");
        std::vector<std::string> n2;
        d.rotors["rotor"].aero->channel_names(n2);
        check(!d.rotors["rotor"].aero->unsteady() && n2.size() + 1==n.size(), "without them: quasi-steady, no channel tau1");
        check(err(head + " dynamic_stall xyz\n").find("dynamic_stall must be oye or hgm") != std::string::npos, "dynamic_stall xyz: an error");
        std::remove("ua_test_af.dat"); std::remove("ua_test_blade.dat");
    }
}

int main()
{
    test_lag();
    test_inflow();
    test_polar();
    test_stall();
    test_rotor();
    std::cout<<(nfail==0 ? "all passed" : std::to_string(nfail) + " failed")<<std::endl;
    return nfail==0 ? 0 : 1;
}
