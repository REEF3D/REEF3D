// REEF3D::CORAL mooring against MoorDyn v2 (no MPI, no REEF3D, no MoorDyn): the cases of
// Regression/moordyn/coral_moordyn_cases.h run in CORAL and are compared with the MoorDyn time series in
// Regression/moordyn/ref/ (written by Regression/moordyn/coral_moordyn_compare --ref); per case the largest
// difference of every channel (mooring force and moment on the body, end tensions of the lines, positions of
// the free points, and for a free body its motion) relative to the channel's range must stay below the case's
// tolerance (snap loads: the mean peaks per period). Published and measured results: the Lowe and Langley (2006)
// jumper (top tension), the chain of Azcona et al. (2017) in the basin of Ecole Centrale de Nantes (static
// fairlead tension, the dynamic peak against the quasi-static one). The cases marked 'unit' run; any case by name. Then the MoorDyn-type input options of the dynamic line: BA < 0 (damping ratio), EA_curve
// (stress-strain table), current_profile, outputs (node positions and segment tensions).
// Build:  g++ -O2 -DEIGEN_NO_DEBUG -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         coral_moordyn_test.cpp <CORAL objects> -o coral_moordyn_test   (or coral_tests.sh)
// Run:    ./coral_moordyn_test [case ...]   from Regression/unit
// Architect: Hans Bihs
#include"../moordyn/coral_moordyn_cases.h"
#include"coral/io/coral_input.h"
#include<algorithm>
#include<cstdio>
#include<filesystem>
#include<fstream>
#include<iostream>
#include<string>

using namespace cmd;

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

static void test_cases(const std::vector<std::string> &only)
{
    for(const case_spec &c : cases())
    {
        if(only.empty() ? !c.unit : std::find(only.begin(), only.end(), c.name)==only.end())
            continue;
        std::cout<<c.name<<": "<<c.title<<std::endl;
        series ref;
        try {ref = series::read("../moordyn/ref/" + c.name + ".csv");}
        catch(const std::exception&) {check(false, "no reference ../moordyn/ref/" + c.name + ".csv"); continue;}
        series s;
        try {s = run_coral(c);}
        catch(const std::exception &e) {check(false, std::string("CORAL: ") + e.what()); continue;}
        double worst = 0.0;
        std::string wn;
        for(const auto &e : c.peaks ? compare_peaks(ref, s, c.t_cmp, c.T) : compare(ref, s, c.t_cmp))
        {
            if(!c.free_body && (e.name=="x" || e.name=="z" || e.name=="pitch"))
                continue;
            if(e.rel > worst) {worst = e.rel; wn = e.name;}
        }
        check(worst <= c.tol, fmt(c.peaks ? "largest difference of the mean peaks to MoorDyn %.2f %%" : "largest difference to MoorDyn %.2f %% of the range",
                                  100.0*worst) + " (" + wn + "), tolerance " + fmt("%.1f %%", 100.0*c.tol));
        // published or measured results
        if(!c.ext_ref.empty())
        {
            try
            {
                const auto ext = read_ext("../moordyn/" + c.ext_ref);
                double md, rms, tm;
                compare_ext(s, ext, c.ext_channel, c.ext_scale, c.ext_t0, md, rms, tm);
                check(md <= c.ext_tol, fmt("against the published solution: largest difference %.3f (at %.1f s), rms %.3f, tolerance %.2f", md, tm, rms, c.ext_tol));
            }
            catch(const std::exception &e) {check(false, e.what());}
        }
        if(std::isfinite(c.ext_static))
        {
            const double v = c.ext_scale*s.rows[0][s.col(c.ext_channel)];
            check(std::fabs(v/c.ext_static - 1.0) <= c.ext_static_tol,
                  fmt("static: measured %.4g, CORAL %.4g (%+.2f %%), tolerance %.1f %%", c.ext_static, v, 100.0*(v/c.ext_static - 1.0), 100.0*c.ext_static_tol));
        }
        if(c.peaks && std::isfinite(c.ext_peak))
        {
            const double v = c.ext_scale*mean_peak(s, c.ext_channel, c.t_cmp, c.T);
            check(std::fabs(v/c.ext_peak - 1.0) <= c.ext_peak_tol,
                  fmt("mean peak: published %.4g, CORAL %.4g (%+.1f %%), tolerance %.0f %%", c.ext_peak, v, 100.0*(v/c.ext_peak - 1.0), 100.0*c.ext_peak_tol));
        }
    }
}

static void test_options()
{
    std::cout<<"MoorDyn-type options of the dynamic line"<<std::endl;
    // BA < 0: the damping ratio of a segment's axial mode, BA = -BA l sqrt(EA m)
    {
        coral::model m;
        coral::mooring_system s(m);
        coral::mooring_props p;
        p.L = 100.0; p.d = 0.1; p.mass = 50.0; p.EA = 1.0e8; p.BA = -0.8;
        coral::mooring_end a, b;
        a.p = Vec3(0.0, 0.0, -50.0);
        b.p = Vec3(99.0, 0.0, -50.0);
        s.add_line(m, "l", a, b, p, 10);
        const double BA = 0.8*10.0*std::sqrt(1.0e8*50.0);
        check(std::fabs(s.line_props(0).BA - BA) < 1.0e-9*BA, fmt("BA -0.8: %.6g N s (0.8 l sqrt(EA m) = %.6g N s)", s.line_props(0).BA, BA));
    }
    // stress-strain table: linear in between, the last slope beyond, the secant and the largest tangent stiffness
    {
        coral::mooring_props p;
        p.EA_strain = {0.0, 0.01, 0.03};
        p.EA_tension = {0.0, 1.0e6, 5.0e6};
        check(std::fabs(p.elastic(0.005) - 0.5e6) < 1e-6 && std::fabs(p.elastic(0.02) - 3.0e6) < 1e-6
              && std::fabs(p.elastic(0.04) - 7.0e6) < 1e-6,
              fmt("EA_curve: T(0.005) %.3g N, T(0.02) %.3g N, T(0.04) %.3g N (beyond the table)", p.elastic(0.005), p.elastic(0.02), p.elastic(0.04)));
        check(std::fabs(p.EA_secant() - 5.0e6/0.03) < 1e-3 && std::fabs(p.EA_max() - 2.0e8) < 1e-3,
              fmt("EA_curve: secant %.4g N, largest tangent %.4g N", p.EA_secant(), p.EA_max()));
    }
    // current profile: linear in z, the end values outside, plus the uniform current
    {
        coral::model m;
        coral::mooring_system s(m, nullptr, Vec3(0.1, 0.0, 0.0));
        s.set_current_profile({-100.0, 0.0}, {Vec3(0.2, 0.0, 0.0), Vec3(1.2, 0.5, 0.0)});
        const Vec3 u0 = s.current_at(-150.0), u1 = s.current_at(-50.0), u2 = s.current_at(10.0);
        check(std::fabs(u0(0) - 0.3) < 1e-12 && std::fabs(u1(0) - 0.8) < 1e-12 && std::fabs(u1(1) - 0.25) < 1e-12
              && std::fabs(u2(0) - 1.3) < 1e-12, fmt("current_profile: u(-150) %.2f, u(-50) %.2f %.3f, u(10) %.2f m/s", u0(0), u1(0), u1(1), u2(0)));
    }
    // the input keywords: EA_curve and current_profile files, BA < 0, outputs p and t
    {
        const std::string dir = (std::filesystem::temp_directory_path()/"coral_moordyn_test").string();
        std::filesystem::create_directories(dir);
        {
            std::ofstream f(dir + "/curve.txt");
            f<<"# strain tension\n0 0\n0.01 1.0e6\n0.03 5.0e6\n";
            std::ofstream g(dir + "/current.txt");
            g<<"# z ux uy uz\n-50 0.2 0 0\n0 1.0 0 0\n";
        }
        const std::string text =
            "body buoy mass 1e5 inertia 1e6 1e6 1e6 cog 0 0 0\n"
            "line_type poly d 0.16 mass 26 EA_curve curve.txt BA -1.0\n"
            "mooring L1 dynamic a ground -150 0 -50 b buoy -5 0 -5 type poly L 150 nseg 10 current_profile current.txt outputs pt\n";
        coral::coral_case cc;
        bool ok = true;
        try {cc.read_string(text, dir);}
        catch(const std::exception &e) {std::cout<<"  "<<e.what()<<std::endl; ok = false;}
        check(ok, "input: EA_curve, BA -1.0, current_profile, outputs pt");
        if(ok)
        {
            const coral::mooring_system *s = nullptr;
            for(const auto &l : cc.m.load_list)
                if(auto q = dynamic_cast<const coral::mooring_system*>(l.get())) s = q;
            check(s != nullptr, "input: the dynamic line");
            if(s)
            {
                const coral::mooring_props &p = s->line_props(0);
                const double BA = 15.0*std::sqrt(p.EA_secant()*26.0);
                check(p.EA_strain.size()==3 && std::fabs(p.BA - BA) < 1e-9*BA && std::fabs(s->current_at(-25.0)(0) - 0.6) < 1e-12,
                      fmt("input: 3 table points, BA %.4g N s (%.4g), current at -25 m %.2f m/s", p.BA, BA, s->current_at(-25.0)(0)));
                std::vector<std::string> ch;
                s->channel_names(ch);
                const bool has_p = std::find(ch.begin(), ch.end(), "L1_n10_z") != ch.end();
                const bool has_t = std::find(ch.begin(), ch.end(), "L1_T9") != ch.end() && std::find(ch.begin(), ch.end(), "L1_T10")==ch.end();
                std::vector<double> v;
                s->channel_values(cc.m, 0.0, v);
                check(has_p && has_t && v.size()==ch.size(), fmt("outputs pt: %g channels (node positions and segment tensions)", (double)ch.size()));
            }
        }
        bool bad = false;
        try {coral::coral_case c2; c2.read_string("body buoy mass 1 inertia 1 1 1\nmooring L1 dynamic a ground -10 0 -5 b buoy 0 0 0 d 0.1 mass 5 EA 1e6 L 11 outputs px\n", dir);}
        catch(const std::exception&) {bad = true;}
        check(bad, "input: outputs with an unknown flag is refused");
    }
}

int main(int argc, char **argv)
{
    std::vector<std::string> only;
    for(int i=1; i<argc; ++i) only.push_back(argv[i]);
    if(only.empty())
        test_options();
    test_cases(only);
    std::cout<<(nfail ? "FAILED: " + std::to_string(nfail) : std::string("all passed"))<<std::endl;
    return nfail ? 1 : 0;
}
