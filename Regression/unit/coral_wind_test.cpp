// REEF3D::CORAL unit test of the wind inflow (coral/wind): steady wind with power-law and log shear, veer,
// upflow, ramp and the hub position; uniform wind in time and InflowWind's uniform wind file; IEC 61400-1 classes,
// turbulence models and the deterministic events EOG, EDC, ECD, EWS; IEC Kaimal turbulence by Veers' method
// (standard deviations, spectra, coherence of u, no coherence of v and w by default, seeds, scaling at the hub);
// TurbSim .bts written and read back, frozen turbulence, tower points; the keywords and the wind probe.
// Architect: Hans Bihs
#include"coral/wind/coral_wind.h"
#include"coral/wind/coral_wind_grid.h"
#include"coral/wind/coral_wind_probe.h"
#include"coral/io/coral_input.h"
#include<cmath>
#include<complex>
#include<cstdio>
#include<fstream>
#include<iostream>
#include<string>
#include<unsupported/Eigen/FFT>

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

static void test_steady()
{
    std::cout<<"steady wind"<<std::endl;
    wind_field w(10.0, 90.0, 30.0*deg);
    w.set_power_law(0.2);
    double e = 0.0;
    for(double z : {10.0, 50.0, 90.0, 150.0})
    {
        const Vec3 v = w.velocity(Vec3(5.0, -3.0, z), 1.0);
        const double U = 10.0*std::pow(z/90.0, 0.2);
        e = std::max(e, (v - Vec3(U*std::cos(30*deg), U*std::sin(30*deg), 0.0)).norm());
    }
    check(e < 1e-12, fmt("power law 0.2, direction 30 deg (err %.1e)", e));
    check(w.velocity(Vec3(0, 0, -1.0), 0.0).norm()==0.0, "no wind below z = 0");
    wind_field l(10.0, 90.0);
    l.set_log_law(0.03);
    const double ul = l.velocity(Vec3(0, 0, 30.0), 0.0)(0);
    check(std::fabs(ul - 10.0*std::log(30.0/0.03)/std::log(90.0/0.03)) < 1e-12, fmt("log law z0 0.03 m: U(30 m) %.4f", ul));
    wind_field v(10.0, 90.0);
    v.set_veer(0.1*deg);
    v.set_upflow(8.0*deg);
    const Vec3 a = v.velocity(Vec3(0, 0, 140.0), 0.0);
    const double th = std::atan2(a(1), a(0)), up = std::atan2(a(2), std::hypot(a(0), a(1)));
    check(std::fabs(th - 5.0*deg) < 1e-12 && std::fabs(up - 8.0*deg) < 1e-12 && std::fabs(a.norm() - 10.0) < 1e-12,
          fmt("veer 0.1 deg/m: 5 deg at 50 m above z_ref (%.4f deg), upflow 8 deg (%.4f), |U| kept", th/deg, up/deg));
    v.set_ramp(10.0);
    check(std::fabs(v.velocity(Vec3(0, 0, 90.0), 5.0).norm() - 5.0) < 1e-12 && v.velocity(Vec3(0, 0, 90.0), 0.0).norm()==0.0,
          "ramp: half at t_ramp/2, zero at 0");
}

static void test_series()
{
    std::cout<<"uniform wind in time"<<std::endl;
    {
        std::ofstream f("wind_test_uniform.hh");
        f<<"! InflowWind uniform wind\n! Time Wind Wind Vert. Horiz. Pwr.Law Lin.Vert. Gust\n";
        f<<"0.0  8.0  0.0  0.0  0.0  0.0  0.0  0.0\n";
        f<<"10.0 12.0 20.0 0.5  0.1  0.2  0.0  1.0\n";
        f<<"20.0 12.0 20.0 0.5  0.0  0.14 0.2  0.0\n";
    }
    std::vector<wind_field::series_row> rows;
    std::string msg;
    const bool ok = wind_field::read_uniform("wind_test_uniform.hh", rows, msg);
    check(ok && rows.size()==3, "read InflowWind uniform wind file " + msg);
    if(!ok) return;
    std::remove("wind_test_uniform.hh");
    wind_field w(0.0, 90.0);
    w.set_series(rows, 126.0);
    // t = 5: U 10, direction -10 deg (clockwise from above in the file), w 0.25, hshear 0.05, alpha 0.1, gust 0.5
    const Vec3 p(0.0, 20.0, 60.0);
    const double eta = 20.0*std::cos(-10*deg);    // lateral position in the wind frame (direction -10 deg)
    const double U = 10.0*(std::pow(60.0/90.0, 0.1) + 0.05*eta/126.0) + 0.5;
    const Vec3 v = w.velocity(p, 5.0);
    const Vec3 ex(U*std::cos(-10*deg), U*std::sin(-10*deg), 0.25);
    check((v - ex).norm() < 1e-9, fmt("t 5 s: speed, direction, vertical speed, shears, gust interpolated (err %.1e)", (v - ex).norm()));
    const double U20 = 12.0*(std::pow(110.0/90.0, 0.14) + 0.2*20.0/126.0);
    check(std::fabs(w.speed(0.0, 110.0, 25.0) - U20) < 1e-12, "after the last row: held; linear vertical shear");
    check(std::fabs(w.hub_speed(10.0) - 13.0) < 1e-12, "hub speed with the gust");
}

static void test_iec()
{
    std::cout<<"IEC 61400-1 classes and events"<<std::endl;
    iec_turbine tc;
    check(tc.set_class("I", "B") && tc.V_ref==50.0 && tc.I_ref==0.14 && !tc.set_class("IV", "A"), "class I B: V_ref 50, I_ref 0.14");
    tc.z_hub = 90.0; tc.D = 126.0;
    const double s = tc.sigma_ntm(11.4);
    check(std::fabs(s - 0.14*(0.75*11.4 + 5.6)) < 1e-12, fmt("NTM sigma_1 at 11.4 m/s %.4f m/s", s));
    check(std::fabs(tc.sigma_etm(11.4) - 2.0*0.14*(0.072*(10.0/2.0 + 3.0)*(11.4/2.0 - 4.0) + 10.0)) < 1e-12, fmt("ETM %.4f m/s", tc.sigma_etm(11.4)));
    check(std::fabs(tc.V_e50(90.0) - 70.0) < 1e-12 && std::fabs(tc.V_e1(90.0) - 56.0) < 1e-12, "V_e50 = 1.4 V_ref, V_e1 = 0.8 V_e50 at the hub");
    check(tc.Lambda1()==42.0, "Lambda_1 42 m above 60 m");

    const double V = 11.4;
    wind_field w(V, 90.0);
    w.set_power_law(0.2);
    // EOG
    w.set_event(iec_event::eog, tc, 20.0);
    const double Vg = std::min(1.35*(56.0 - V), 3.3*s/(1.0 + 0.1*126.0/42.0));
    double e = 0.0, vmin = 1e9, vmax = -1e9;
    for(double t=20.0; t<=30.5; t+=0.05)
    {
        const double tau = t - 20.0;
        const double ex = V - 0.37*Vg*std::sin(3*M_PI*tau/10.5)*(1 - std::cos(2*M_PI*tau/10.5));
        const double u = w.velocity(Vec3(0, 0, 90.0), t)(0);
        e = std::max(e, std::fabs(u - ex));
        vmin = std::min(vmin, u); vmax = std::max(vmax, u);
    }
    check(e < 1e-9, fmt("EOG: V_gust %.3f m/s, formula of 6.3.3.2 (err %.1e), dip to %.2f, peak %.2f m/s", Vg, e, vmin, vmax));
    check(std::fabs(w.velocity(Vec3(0, 0, 90.0), 19.9)(0) - V) < 1e-12 && std::fabs(w.velocity(Vec3(0, 0, 90.0), 31.0)(0) - V) < 1e-12,
          "EOG: steady before and after");
    // EDC
    w.set_event(iec_event::edc_pos, tc, 0.0);
    const double the = 4.0*std::atan(s/(V*(1.0 + 0.1*126.0/42.0)));
    const Vec3 a = w.velocity(Vec3(0, 0, 90.0), 3.0), b = w.velocity(Vec3(0, 0, 90.0), 7.0);
    check(std::fabs(std::atan2(a(1), a(0)) - 0.5*the) < 1e-12 && std::fabs(std::atan2(b(1), b(0)) - the) < 1e-12,
          fmt("EDC+: half of theta_e %.2f deg after 3 s, all of it after 6 s", the/deg));
    w.set_event(iec_event::edc_neg, tc, 0.0);
    check(std::fabs(std::atan2(w.velocity(Vec3(0, 0, 90.0), 7.0)(1), w.velocity(Vec3(0, 0, 90.0), 7.0)(0)) + the) < 1e-12, "EDC-: the other way");
    // ECD
    w.set_event(iec_event::ecd_pos, tc, 0.0);
    const Vec3 c = w.velocity(Vec3(0, 0, 90.0), 12.0);
    check(std::fabs(c.norm() - (V + 15.0)) < 1e-12 && std::fabs(std::atan2(c(1), c(0)) - 720.0/V*deg) < 1e-12,
          fmt("ECD: +15 m/s and 720/V_hub = %.1f deg after 10 s", 720.0/V));
    // EWS
    w.set_event(iec_event::ews_v_pos, tc, 0.0);
    const double dU = 0.5*(2.5 + 0.2*6.4*s*std::pow(126.0/42.0, 0.25))*2.0;
    const double ut = w.velocity(Vec3(0, 0, 153.0), 6.0)(0), ub = w.velocity(Vec3(0, 0, 27.0), 6.0)(0);
    check(std::fabs(ut - (V*std::pow(1.7, 0.2) + dU)) < 1e-12 && std::fabs(ub - (V*std::pow(0.3, 0.2) - dU)) < 1e-12,
          fmt("EWS vertical: +-%.3f m/s at the rotor's top and bottom at T/2", dU));
    w.set_event(iec_event::ews_h_neg, tc, 0.0);
    check(std::fabs(w.velocity(Vec3(0, 63.0, 90.0), 6.0)(0) - (V - dU)) < 1e-12, "EWS horizontal-: -dU at eta = +D/2");
    check(iec_event_from("ews_h-")==iec_event::ews_h_neg && iec_event_from("x")==iec_event::none, "event names");
}

// statistics of a fluctuation grid
static void test_kaimal()
{
    std::cout<<"IEC Kaimal turbulence (Veers)"<<std::endl;
    iec_turbine tc;
    tc.set_class("I", "B"); tc.z_hub = 90.0; tc.D = 126.0;
    const double V = 11.4, s1 = tc.sigma_ntm(V);
    wind_grid::kaimal_settings ks;
    ks.ny = 7; ks.nz = 7; ks.width = 120.0; ks.height = 120.0; ks.T = 3200.0; ks.dt = 0.2; ks.seed = 3;
    const wind_grid g = wind_grid::kaimal(tc, V, s1, ks);
    check(g.nt==16384 && g.periodic && !g.has_mean && std::fabs(g.dy - 20.0) < 1e-12 && std::fabs(g.z_bottom - 30.0) < 1e-12,
          fmt("grid 7 x 7, 20 m apart from z = 30 m, %.0f levels of 0.2 s (2^14, periodic)", g.nt));
    const int P = g.ny*g.nz, N = g.nt;
    auto x = [&](int it, int j, int c) {return (double)g.data[((size_t)it*P + j)*3 + c];};
    // standard deviations over all points
    double var[3] = {0, 0, 0}, mean[3] = {0, 0, 0};
    for(int c=0; c<3; ++c)
    {
        for(int j=0; j<P; ++j)
            for(int it=0; it<N; ++it) {mean[c] += x(it, j, c); var[c] += x(it, j, c)*x(it, j, c);}
        mean[c] /= (double)N*P; var[c] /= (double)N*P;
    }
    // expected variance of the series: the Kaimal spectrum summed over the lines of the series
    double ev[3] = {0, 0, 0};
    const double df = 1.0/(N*g.dt);
    for(int c=0; c<3; ++c)
        for(int k=1; k<N/2; ++k)
            ev[c] += wind_grid::kaimal_spectrum(tc, c, k*df, V, s1)*df;
    check(std::fabs(mean[0]) < 1e-3 && std::fabs(std::sqrt(var[0]/ev[0]) - 1.0) < 0.05 && std::fabs(std::sqrt(var[1]/ev[1]) - 1.0) < 0.05
          && std::fabs(std::sqrt(var[2]/ev[2]) - 1.0) < 0.05,
          fmt("sigma u v w %.3f %.3f %.3f m/s, of the spectra %.3f", std::sqrt(var[0]), std::sqrt(var[1]), std::sqrt(var[2]), std::sqrt(ev[0])));
    check(std::fabs(std::sqrt(ev[0]) - s1) < 0.03*s1 && std::fabs(std::sqrt(ev[2]) - 0.5*s1) < 0.03*s1,
          fmt("the spectra's sigma_1 %.3f and sigma_3 %.3f of sigma_1 %.3f (lines 1/T .. f_Nyquist)", std::sqrt(ev[0]), std::sqrt(ev[2]), s1));
    // spectra and cross-spectra by the FFT of the whole series
    Eigen::FFT<double> fft;
    std::vector<std::vector<std::vector<std::complex<double>>>> F(3, std::vector<std::vector<std::complex<double>>>(P));
    for(int c=0; c<3; ++c)
        for(int j=0; j<P; ++j)
        {
            std::vector<double> s(N);
            for(int it=0; it<N; ++it) s[it] = x(it, j, c);
            fft.fwd(F[c][j], s);
        }
    auto band = [&](int c, int a, int b, double f0, double f1, double &S, double &Sx)
    {
        S = Sx = 0.0;
        for(int k=(int)(f0/df); k<=(int)(f1/df); ++k)
        {
            S += 0.5*(std::norm(F[c][a][k]) + std::norm(F[c][b][k]));
            Sx += std::real(F[c][a][k]*std::conj(F[c][b][k]));
        }
    };
    // u spectrum in octave bands against the Kaimal spectrum (all points)
    double emax = 0.0;
    for(double f0 : {0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64})
    {
        double S = 0.0, Se = 0.0;
        for(int j=0; j<P; ++j)
            for(int k=(int)(f0/df); k<(int)(2*f0/df); ++k)
                S += std::norm(F[0][j][k])*2.0/((double)N*N)/df/P;
        for(int k=(int)(f0/df); k<(int)(2*f0/df); ++k)
            Se += wind_grid::kaimal_spectrum(tc, 0, k*df, V, s1);
        emax = std::max(emax, std::fabs(S/Se - 1.0));
    }
    check(emax < 0.1, fmt("u spectrum in octave bands 0.01 - 1.28 Hz = Kaimal within %.1f %%", 100*emax));
    // coherence of u between neighbours (20 m) and across the grid (120 m), v and w uncorrelated
    auto coh = [&](int c, int a, int b, double f0, double f1)
    {
        double S, Sx, St = 0, Sxt = 0;
        band(c, a, b, f0, f1, S, Sx);
        St += S; Sxt += Sx;
        return Sxt/St;
    };
    double ce = 0.0;
    std::string cs;
    for(double f0 : {0.005, 0.01, 0.02, 0.04})
    {
        // all neighbour pairs along y and z, averaged; IEC's coherence averaged over the lines of the band
        double c20 = 0.0, cth = 0.0;
        int n = 0;
        for(int i=0; i<g.nz; ++i)
            for(int j=0; j<g.ny; ++j)
            {
                if(j+1 < g.ny) {c20 += coh(0, i*g.ny + j, i*g.ny + j + 1, f0, 2*f0); ++n;}
                if(i+1 < g.nz) {c20 += coh(0, i*g.ny + j, (i+1)*g.ny + j, f0, 2*f0); ++n;}
            }
        int nk = 0;
        for(int k=(int)(f0/df); k<=(int)(2*f0/df); ++k, ++nk)
            cth += wind_grid::iec_coherence(tc, 20.0, k*df, V);
        ce = std::max(ce, std::fabs(c20/n - cth/nk));
        cs += fmt(" %.3f/%.3f", c20/n, cth/nk);
    }
    check(ce < 0.04, "u coherence of neighbours 20 m apart, octaves from 0.005 Hz (grid/IEC):" + cs);
    double cv = 0.0, cw = 0.0;
    int nv = 0;
    for(int i=0; i<g.nz; ++i)
        for(int j=0; j+1<g.ny; ++j, ++nv)
        {
            cv += coh(1, i*g.ny + j, i*g.ny + j + 1, 0.005, 0.05);
            cw += coh(2, i*g.ny + j, i*g.ny + j + 1, 0.005, 0.05);
        }
    check(std::fabs(cv/nv) < 0.03 && std::fabs(cw/nv) < 0.03, fmt("v and w without spatial coherence by default (mean %.3f, %.3f)", cv/nv, cw/nv));
    // seeds, coherence of v and w, scaling at the hub
    wind_grid::kaimal_settings k2 = ks;
    k2.T = 400.0;
    const wind_grid g1 = wind_grid::kaimal(tc, V, s1, k2), g2 = wind_grid::kaimal(tc, V, s1, k2);
    k2.seed = 4;
    const wind_grid g3 = wind_grid::kaimal(tc, V, s1, k2);
    check(g1.data==g2.data && g1.data != g3.data, "the same seed gives the same field, another seed another");
    k2.seed = 3; k2.scale_hub = true; k2.coherence_vw = true;
    const wind_grid g4 = wind_grid::kaimal(tc, V, s1, k2);
    const int jh = 3*7 + 3;
    double sd[3] = {0, 0, 0};
    for(int c=0; c<3; ++c)
    {
        double m = 0, q = 0;
        for(int it=0; it<g4.nt; ++it) {const double y = g4.data[((size_t)it*P + jh)*3 + c]; m += y; q += y*y;}
        m /= g4.nt;
        sd[c] = std::sqrt(q/g4.nt - m*m);
    }
    check(std::fabs(sd[0] - s1) < 1e-4 && std::fabs(sd[1] - 0.8*s1) < 1e-4 && std::fabs(sd[2] - 0.5*s1) < 1e-4,
          fmt("scale_hub: sigma at the hub %.4f %.4f %.4f exactly", sd[0], sd[1], sd[2]));
}

static void test_bts()
{
    std::cout<<"TurbSim .bts, frozen turbulence, tower points"<<std::endl;
    iec_turbine tc;
    tc.set_class("II", "A"); tc.z_hub = 90.0; tc.D = 126.0;
    wind_grid::kaimal_settings ks;
    ks.ny = 5; ks.nz = 6; ks.width = 100.0; ks.height = 100.0; ks.T = 100.0; ks.dt = 0.25; ks.seed = 9;
    auto g = std::make_shared<wind_grid>(wind_grid::kaimal(tc, 10.0, tc.sigma_ntm(10.0), ks));
    wind_field w(10.0, 90.0);
    w.set_power_law(0.14);
    std::string msg;
    check(g->write_bts("wind_test.bts", "CORAL unit test", msg, &w), "write " + msg);
    wind_grid r;
    check(r.read_bts("wind_test.bts", msg), "read back " + msg);
    std::remove("wind_test.bts");
    check(r.ny==5 && r.nz==6 && r.nt==g->nt && r.periodic && r.has_mean && std::fabs(r.u_hub - 10.0) < 1e-6
          && std::fabs(r.z_bottom - 40.0) < 1e-5 && std::fabs(r.dy - 25.0) < 1e-5 && std::fabs(r.dz - 20.0) < 1e-5 && std::fabs(r.dt - 0.25) < 1e-7,
          "header: 5 x 6 points 25 / 20 m apart from z 40 m, levels, dt, hub speed, periodic (identifier 8)");
    // the values: mean + fluctuation within the 16-bit resolution
    double e = 0.0, range = 0.0;
    for(int it=0; it<r.nt; it+=7)
        for(int i=0; i<r.nz; ++i)
            for(int j=0; j<r.ny; ++j)
            {
                const size_t q = (((size_t)it*r.nz + i)*r.ny + j)*3;
                const double U = w.speed(r.eta0() + j*r.dy, r.z_bottom + i*r.dz, 0.0);
                e = std::max(e, (double)std::fabs(r.data[q] - (g->data[q] + U)));
                e = std::max(e, (double)std::fabs(r.data[q+2] - g->data[q+2]));
                range = std::max(range, (double)std::fabs(r.data[q]));
            }
    check(e < 2.0*range/65535.0 + 1e-6, fmt("velocities (with the shear profile) to the 16-bit resolution (err %.1e m/s)", e));
    // frozen turbulence: a point dx downwind sees the field dx/u_hub later
    r.origin = Vec3(100.0, 50.0, 0.0);
    r.direction = 0.0;
    const Vec3 a = r.velocity(Vec3(100.0, 50.0, 95.0), 10.0), b = r.velocity(Vec3(120.0, 50.0, 95.0), 12.0);
    check((a - b).norm() < 1e-6, "frozen: 20 m downwind 2 s later the same velocity");
    r.direction = 90.0*deg;
    const Vec3 c = r.velocity(Vec3(100.0, 70.0, 95.0), 12.0);
    check(std::fabs(c(1) - a(0)) < 1e-5 && std::fabs(c(0) + a(1)) < 1e-5, "direction 90 deg: the field turned with it");
    // periodic in time
    r.direction = 0.0;
    const Vec3 d = r.velocity(Vec3(100.0, 50.0, 95.0), 10.0 + r.period());
    check((d - a).norm() < 1e-5, fmt("periodic after %.0f s", r.period()));
    // tower points: a 1 x 1 grid at z 50 with two tower points
    wind_grid t;
    t.ny = 1; t.nz = 1; t.nt = 1; t.ntower = 2; t.dz = 10.0; t.z_bottom = 50.0; t.u_hub = 5.0;
    t.data = {6.0f, 0.0f, 0.0f};
    t.tower = {4.0f, 0.0f, 0.0f, 3.0f, 0.0f, 0.0f};
    check(std::fabs(t.local(0.0, 45.0, 0.0)(0) - 5.0) < 1e-6 && std::fabs(t.local(0.0, 35.0, 0.0)(0) - 3.5) < 1e-6
          && std::fabs(t.local(0.0, 15.0, 0.0)(0) - 1.5) < 1e-6 && std::fabs(t.local(30.0, 60.0, 0.0)(0) - 6.0) < 1e-6,
          "below the grid: tower points 10 m apart, then to zero at z = 0; outside: the edge");
    // a fluctuation grid on a mean wind
    wind_field m(10.0, 90.0);
    m.set_turbulence(g);
    const Vec3 f = m.velocity(Vec3(0.0, 10.0, 80.0), 7.0);
    const Vec3 fl = g->local(10.0, 80.0, 7.0);
    check((f - (Vec3(10.0, 0.0, 0.0) + fl)).norm() < 1e-9, "wind_field + turbulence = mean + the grid's fluctuation");
}

static void test_input()
{
    std::cout<<"keywords"<<std::endl;
    coral_case c;
    c.read_string("time dt 0.1 end 1\n"
                  "body nacelle mass 1000 cog 0 0 90\n"
                  "joint fix fixed ground nacelle point 0 0 90\n"
                  "rho_air 1.2\n"
                  "wind steady U 11.4 shear 0.2 direction 10\n"
                  "wind iec I B z_hub 90 D 126\n"
                  "wind event eog t0 0.2\n"
                  "wind turbulence kaimal ny 5 nz 5 T 50 dt 0.1 seed 2\n"
                  "wind_probe hub 0 0 90 body nacelle\n", ".");
    check(c.wind != nullptr && c.rho_air==1.2, "wind steady + iec + event + turbulence, rho_air");
    auto *w = dynamic_cast<wind_field*>(c.wind.get());
    check(w && w->z_ref()==90.0 && w->iec().I_ref==0.14, "z_ref from z_hub of wind iec");
    load *probe = nullptr;
    for(auto &l : c.m.load_list)
        if(l->name()=="hub") probe = l.get();
    std::vector<std::string> n;
    std::vector<double> v;
    if(probe) {probe->channel_names(n); probe->channel_values(c.m, 0.5, v);}
    const Vec3 u = c.wind->velocity(Vec3(0, 0, 90.0), 0.5);
    check(probe && n.size()==4 && n[0]=="hub_u" && v.size()==4 && std::fabs(v[0] - u(0)) < 1e-12 && std::fabs(v[3] - u.norm()) < 1e-12,
          fmt("wind_probe channels hub_u .. hub_U (U %.3f m/s)", u.norm()));
    auto err = [](const std::string &s)
    {
        try {coral_case k; k.read_string(s, ".");} catch(const std::exception &e) {return std::string(e.what());}
        return std::string();
    };
    check(err("wind event eog t0 1\nwind steady U 10 z_ref 90\n").find("needs wind iec") != std::string::npos, "event without iec: error");
    check(err("body b mass 1\nwind_probe p 0 0 1\n").find("needs a wind") != std::string::npos, "probe without wind: error");
    check(err("wind steady U 10\n").find("z_ref") != std::string::npos, "no reference height: error");
    check(err("wind turbsim nofile.bts\n").find("cannot open") != std::string::npos, "missing TurbSim file: error");
}

int main()
{
    test_steady();
    test_series();
    test_iec();
    test_kaimal();
    test_bts();
    test_input();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
