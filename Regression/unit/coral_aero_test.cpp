// REEF3D::CORAL unit test of the aerodynamics (coral/aero): airfoil polars (AeroDyn's airfoil file with two tables
// and the unsteady parameters, plain tables, interpolation in alpha and log(Re), wrap-around of a -180 .. 180 table);
// blade element momentum (Ning's method): the momentum balance of thrust and torque at the converged element, Buhl's
// high-induction thrust, 2000 random elements converged with a zero residual and the inflow angle relation, no load at the tip, the skew correction; the Glauert
// optimum rotor (BEM's induction = the optimum, Cp = Glauert's Cp_max at the tip speed ratio); rotor thrust, torque,
// power; AeroDyn's blade file.
// Architect: Hans Bihs
#include"coral/aero/coral_airfoil.h"
#include"coral/aero/coral_bem.h"
#include<cmath>
#include<cstdio>
#include<iostream>
#include<random>
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

// linear lift 2 pi alpha up to +-20 deg, constant drag
static airfoil thin(double cd=0.0)
{
    airfoil af;
    af.tables.emplace_back();
    polar_table &p = af.tables.back();
    for(int i=-40; i<=40; ++i)
    {
        const double a = 0.5*i*deg;
        p.alpha.push_back(a); p.cl.push_back(2.0*M_PI*a); p.cd.push_back(cd); p.cm.push_back(0.0);
    }
    return af;
}

static void test_airfoil()
{
    std::cout<<"airfoil polars"<<std::endl;
    const std::string ad =
        "! ------------ AirfoilInfo v1.01.x Input File ----------------------------------\n"
        "! test airfoil\n"
        "\"DEFAULT\"     InterpOrd         ! Interpolation order\n"
        "          1   NonDimArea        ! The non-dimensional area\n"
        "          0   NumCoords         ! coordinates\n"
        "          2   NumTabs           ! Number of airfoil tables in this file.\n"
        "! data for table 1\n"
        "       3.0    Re                ! Reynolds number in millions\n"
        "          0   UserProp\n"
        "True          InclUAdata\n"
        "       -3.2   alpha0            ! 0-lift angle of attack\n"
        "    7.4888    C_nalpha\n"
        "    Default   T_f0\n"
        "          4   NumAlf            ! Number of data lines in the following table\n"
        "!    Alpha      Cl      Cd        Cm\n"
        "   -180.00    0.000   0.0602    0.0000\n"
        "     -10.00  -0.700   0.0200   -0.0100\n"
        "      10.00   1.300   0.0300   -0.0500\n"
        "     180.00    0.000   0.0602    0.0000\n"
        "! data for table 2\n"
        "       0.75   Re\n"
        "          0   UserProp\n"
        "False         InclUAdata\n"
        "          3   NumAlf\n"
        "     -10.00  -0.500   0.0300   -0.0100\n"
        "      10.00   1.100   0.0400   -0.0500\n"
        "      20.00   1.000   0.1000   -0.0800\n";
    airfoil af;
    std::string msg;
    const bool ok = airfoil::read_string(ad, af, msg);
    check(ok && af.tables.size()==2 && af.tables[0].re==0.75e6 && af.tables[1].re==3.0e6, "AeroDyn airfoil file: two tables sorted by Re " + msg);
    if(!ok) return;
    check(af.tables[1].params.at("alpha0")==-3.2 && af.tables[1].params.count("T_f0")==0 && af.tables[1].alpha.size()==4,
          "unsteady parameters kept by name (Default left out), NumAlf rows");
    const airfoil_coef c = af.table_coefficients(1, 5.0*deg);
    check(std::fabs(c.cl - 0.8) < 1e-12 && std::fabs(c.cd - 0.0275) < 1e-12 && std::fabs(c.cm + 0.04) < 1e-12, "linear in alpha");
    const airfoil_coef w = af.table_coefficients(1, 200.0*deg), w2 = af.table_coefficients(1, -160.0*deg);
    check(std::fabs(w.cl - w2.cl) < 1e-12 && std::fabs(w.cd - w2.cd) < 1e-12, "a -180 .. 180 table wraps around");
    const airfoil_coef h = af.table_coefficients(0, 40.0*deg);
    check(std::fabs(h.cl - 1.0) < 1e-12, "a table short of 180 deg is held at its end");
    const double re = std::sqrt(0.75e6*3.0e6);
    const airfoil_coef m = af.coefficients(5.0*deg, re);
    check(std::fabs(m.cl - 0.5*(0.7 + 0.8)) < 1e-12 && std::fabs(af.coefficients(5.0*deg, 1e9).cl - 0.8) < 1e-12,
          "between the tables linear in log(Re), held outside");
    airfoil pl;
    check(airfoil::read_string("# alpha cl cd\n-5 -0.4 0.01\n0 0.2 0.008\n5 0.8 0.012\n", pl, msg) && pl.tables.size()==1
          && std::fabs(pl.coefficients(2.5*deg).cl - 0.5) < 1e-12 && pl.tables[0].cm.empty(), "plain table without cm");
    check(!airfoil::read_string("5 0.8 0.01\n0 0.2 0.01\n", pl, msg) && msg.find("increase") != std::string::npos, "decreasing alpha: error");
}

static void test_element()
{
    std::cout<<"blade element momentum (Ning)"<<std::endl;
    const airfoil af = thin(0.01);
    bem_element e;
    e.r = 20.0; e.R = 40.0; e.R_hub = 2.0; e.B = 3; e.chord = 1.5; e.theta = 4.0*deg; e.af = &af;
    const double rho = 1.225, Vx = 10.0, Vy = 40.0;
    const bem_result r = bem_solve(e, Vx, Vy, rho);
    const double Tm = 4.0*M_PI*e.r*rho*Vx*Vx*r.a*(1.0 - r.a)*r.F, Qm = 4.0*M_PI*e.r*rho*Vx*Vy*r.ap*(1.0 - r.a)*r.F;
    check(r.converged && r.region==1 && r.a > 0.0 && r.a < 1.0/3.0, fmt("converged in the windmill region: a %.4f, a' %.5f, F %.4f", r.a, r.ap, r.F));
    check(std::fabs(e.B*r.fn/Tm - 1.0) < 1e-9 && std::fabs(e.B*r.ft/Qm - 1.0) < 1e-9,
          fmt("momentum balance: B fn = 4 pi r rho Vx^2 a (1-a) F (%.2e), B ft = 4 pi r rho Vx Vy a' (1-a) F (%.2e)", e.B*r.fn/Tm - 1, e.B*r.ft/Qm - 1));
    check(std::fabs(std::tan(r.phi) - Vx*(1 - r.a)/(Vy*(1 + r.ap))) < 1e-10 && std::fabs(r.alpha - (r.phi - e.theta)) < 1e-15,
          fmt("inflow angle %.3f deg, angle of attack %.3f deg", r.phi/deg, r.alpha/deg));
    // heavily loaded: Buhl
    bem_element h = e;
    h.chord = 6.0; h.theta = 1.0*deg;
    const bem_result b = bem_solve(h, Vx, Vy, rho);
    const double CT = h.B*b.fn/(0.5*rho*Vx*Vx*2.0*M_PI*h.r);
    const double CTb = 8.0/9.0 + (4.0*b.F - 40.0/9.0)*b.a + (50.0/9.0 - 4.0*b.F)*b.a*b.a;
    check(b.a > 0.4 && std::fabs(CT - CTb) < 1e-8, fmt("high induction a %.3f: local C_T %.4f = Buhl's %.4f", b.a, CT, CTb));
    // random elements: tan(phi) relation and a zero residual wherever it converged; regions
    std::mt19937 gen(5);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    int n1 = 0, n2 = 0, n3 = 0, bad = 0;
    for(int i=0; i<2000; ++i)
    {
        bem_element q = e;
        q.r = 2.5 + 37.0*u(gen);
        q.chord = 0.3 + 6.0*u(gen);
        q.theta = (-25.0 + 50.0*u(gen))*deg;
        const double vy = 1.0 + 80.0*u(gen);
        const bem_result s = bem_solve(q, Vx, vy, rho);
        if(!s.converged) {++bad; continue;}
        if(s.region==1) ++n1; else if(s.region==2) ++n2; else if(s.region==3) ++n3;
        const double rr = std::fabs(bem_residual(q, Vx, vy, s.phi, bem_options()));
        double tr = std::fabs(std::sin(s.phi)*vy*(1 + s.ap) - std::cos(s.phi)*Vx*(1 - s.a));
        // the tan(phi) relation holds except in the propeller brake with k <= 1 (a = 0 there, Ning)
        if(s.region==2 && s.a==0.0) tr = 0.0;
        if(rr > 1e-7 || tr > 1e-6*(1 + std::fabs(vy*(1 + s.ap)) + std::fabs(Vx*(1 - s.a))))
            ++bad;
    }
    check(bad==0 && n1 > 0 && n2 > 0, fmt("2000 random elements: all converged with zero residual (%.0f windmill, %.0f propeller brake, %.0f beyond pi/2)", n1, n2, n3));
    // tip and skew
    bem_element t = e;
    t.r = t.R;
    const bem_result tr = bem_solve(t, Vx, Vy, rho);
    check(tr.fn==0.0 && tr.ft==0.0 && tr.F==0.0, "no load at the tip with the tip loss");
    const bem_result sk = bem_solve(e, Vx, Vy, rho, bem_options(), 0.2);
    check(std::fabs(sk.a - 1.2*r.a) < 1e-12 && std::fabs(std::tan(sk.phi) - Vx*(1 - sk.a)/(Vy*(1 + sk.ap))) < 1e-12, "skewed wake: a (1 + skew), phi from it");
    bem_options o;
    o.induction = false;
    const bem_result ni = bem_solve(e, Vx, Vy, rho, o);
    check(ni.a==0.0 && std::fabs(ni.W - std::hypot(Vx, Vy)) < 1e-12, "induction off: W = |V|");
}

// Glauert's optimum rotor with wake rotation: lambda_r^2 = (1 - a)(4a - 1)^2/(1 - 3a), a' = (1 - 3a)/(4a - 1)
static double glauert_a(double lr)
{
    double lo = 0.25 + 1e-14, hi = 1.0/3.0 - 1e-14;
    for(int i=0; i<200; ++i)
    {
        const double a = 0.5*(lo + hi);
        const double f = (1 - a)*(4*a - 1)*(4*a - 1)/(1 - 3*a) - lr*lr;
        (f > 0 ? hi : lo) = a;
    }
    return 0.5*(lo + hi);
}

static void test_optimum_rotor()
{
    std::cout<<"Glauert optimum rotor"<<std::endl;
    const airfoil af = thin(0.0);
    const double lam = 7.0, R = 40.0, U = 10.0, Om = lam*U/R, ad = 6.0*deg, cld = 2.0*M_PI*ad;
    rotor_def rd;
    rd.B = 3; rd.R_tip = R; rd.R_hub = 0.0;
    rd.airfoils.push_back(af);
    const int n = 400;
    for(int i=0; i<=n; ++i)
    {
        const double r = R*(0.01 + 0.99*i/n);
        const double lr = lam*r/R;
        const double phi = 2.0/3.0*std::atan(1.0/lr);
        blade_station s;
        s.r = r;
        s.chord = 8.0*M_PI*r*(1.0 - std::cos(phi))/(rd.B*cld);
        s.twist = phi - ad;
        rd.stations.push_back(s);
    }
    bem_options o;
    o.tip_loss = o.hub_loss = false;
    o.drag_in_induction = false;
    const rotor_result res = rotor_steady(rd, U, Om, 0.0, 1.225, o);
    double ea = 0.0, eap = 0.0;
    for(int i=0; i<=n; ++i)
    {
        const double lr = lam*rd.stations[i].r/R;
        const double a = glauert_a(lr);
        ea = std::max(ea, std::fabs(res.el[i].a - a));
        eap = std::max(eap, std::fabs(res.el[i].ap - (1 - 3*a)/(4*a - 1))/((1 - 3*a)/(4*a - 1)));
    }
    check(ea < 1e-6 && eap < 1e-5, fmt("BEM induction = the optimum: a within %.1e, a' within %.1e (relative)", ea, eap));
    // Glauert's Cp_max: 8/lambda^2 int lambda_r^3 a' (1 - a) dlambda_r over the same stations
    double cp = 0.0;
    for(int i=0; i<n; ++i)
    {
        auto g = [&](int j) {const double lr = lam*rd.stations[j].r/R, a = glauert_a(lr); return lr*lr*lr*(1 - 3*a)/(4*a - 1)*(1 - a);};
        cp += 0.5*(g(i) + g(i+1))*lam*(rd.stations[i+1].r - rd.stations[i].r)/R;
    }
    cp *= 8.0/(lam*lam);
    check(std::fabs(res.Cp/cp - 1.0) < 1e-5 && cp > 0.575 && cp < 0.585,
          fmt("Cp %.5f = Glauert's Cp_max %.5f at tip speed ratio 7 (lambda %.2f)", res.Cp, cp, res.lambda));
    check(std::fabs(res.P - res.Q*Om) < 1e-9*res.P && res.Ct > 0.8 && res.Ct < 1.0, fmt("P = Q Omega, C_T %.4f", res.Ct));
    o.tip_loss = true;
    const rotor_result tl = rotor_steady(rd, U, Om, 0.0, 1.225, o);
    check(tl.Cp < res.Cp && tl.Cp > 0.9*res.Cp, fmt("with the tip loss Cp %.4f", tl.Cp));
}

static void test_blade_file()
{
    std::cout<<"AeroDyn blade file"<<std::endl;
    const std::string bl =
        "------- AERODYN v15.00.* BLADE DEFINITION INPUT FILE -------------------------------------\n"
        "test blade\n"
        "======  Blade Properties =================================================================\n"
        "          3   NumBlNds          - Number of blade nodes used in the analysis (-)\n"
        "    BlSpn        BlCrvAC        BlSwpAC        BlCrvAng       BlTwist        BlChord          BlAFID\n"
        "     (m)           (m)            (m)            (deg)         (deg)           (m)              (-)\n"
        "0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  1.3308000E+01  3.5420000E+00  1\n"
        "3.0000000E+01  0.0000000E+00  0.0000000E+00  0.0000000E+00  6.5440000E+00  3.7480000E+00  8\n"
        "6.1500000E+01  0.0000000E+00  0.0000000E+00  0.0000000E+00  1.0600000E-01  1.4190000E+00  8\n";
    std::vector<blade_station> st;
    std::string msg;
    const bool ok = read_aerodyn_blade_string(bl, 1.5, st, msg);
    check(ok && st.size()==3 && std::fabs(st[2].r - 63.0) < 1e-12 && std::fabs(st[1].twist - 6.544*deg) < 1e-12 && st[1].airfoil==7
          && std::fabs(st[0].chord - 3.542) < 1e-12, "3 nodes, r = R_hub + BlSpn, twist, chord, airfoil index " + msg);
    check(!read_aerodyn_blade_string("  4 NumBlNds\n0 0 0 0 1 1 1\n", 1.5, st, msg), "too few nodes: error");
}

int main()
{
    test_airfoil();
    test_element();
    test_optimum_rotor();
    test_blade_file();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
