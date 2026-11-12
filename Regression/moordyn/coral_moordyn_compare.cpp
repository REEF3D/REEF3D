/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/
// REEF3D::CORAL against MoorDyn v2 (MoorDyn-C, linked as a library): every case of coral_moordyn_cases.h in both
// codes, the time series out/<case>_coral.csv and out/<case>_moordyn.csv, the largest difference per channel;
// with --ref the MoorDyn series become the references ref/<case>.csv of the unit test.
//
// MoorDyn side: a coupled body at the body origin carries the body ends (points of type Body1), so the forces
// it returns include the inertia of the end nodes like CORAL's attachments; the anchors are fixed points, the
// free points free. Prescribed motion: position and velocity at the start of each coupling step (MoorDyn moves
// the coupled body linearly over the step). A free body: the planar equations of motion of the case (surge,
// heave, pitch) with the semi-implicit Euler scheme and MoorDyn's force at the end of each coupling step (first
// order; the case's coupling step is small against the body's periods). Same sub-step as CORAL, RK4, the
// stationary initial condition solver.
//
// Build: see README.md (build_moordyn.sh, then make -f Makefile.compare)
// Run:   ./coral_moordyn_compare [--ref] [case ...]
// Architect: Hans Bihs

#include"coral_moordyn_cases.h"
#include"MoorDyn2.h"
#include<cstdio>
#include<filesystem>
#include<iomanip>
#include<iostream>
#include<map>

using namespace cmd;
namespace fs = std::filesystem;

// ------------------------------------------------------------------ MoorDyn input file

static void write_moordyn_input(const case_spec &c, const std::string &dir, double dtM)
{
    fs::create_directories(dir);
    std::ofstream f(dir + "/moordyn.txt");
    f<<std::setprecision(12);
    f<<"--------------------- MoorDyn Input File ------------------------------------\n";
    f<<"CORAL comparison case "<<c.name<<"\n";
    f<<"----------------------- LINE TYPES ------------------------------------------\n";
    f<<"TypeName Diam Mass/m EA BA/-zeta EI Cd Ca CdAx CaAx\n";
    f<<"(name) (m) (kg/m) (N) (N-s/-) (N-m^2) (-) (-) (-) (-)\n";
    for(const auto &t : c.types)
    {
        std::string EA;
        if(t.EA_curve.size() >= 2)
        {
            EA = t.name + "_EA.txt";
            std::ofstream e(dir + "/" + EA);
            e<<std::setprecision(12)<<"strain tension curve of "<<t.name<<"\n(-) (N)\n---\n";
            for(size_t j=0; j<t.EA_curve.size(); ++j)
                e<<t.EA_curve[j][0]<<" "<<t.EA_curve[j][1]<<(j+1 < t.EA_curve.size() ? "\n" : "");
        }
        else
        {
            std::ostringstream s;
            s<<std::setprecision(12)<<t.EA;
            EA = s.str();
        }
        f<<t.name<<" "<<t.d<<" "<<t.mass<<" "<<EA<<" "<<t.BA<<" "<<t.EI<<" "<<t.Cdn<<" "<<t.Can<<" "<<t.Cdt<<" "<<t.Cat<<"\n";
    }
    f<<"---------------------------- BODIES -----------------------------------------------------\n";
    f<<"ID Attachment X0 Y0 Z0 r0 p0 y0 Mass CG* I* Volume CdA* Ca\n";
    f<<"(#) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m) (kg-m^2) (m^3) (m^2) (-)\n";
    const Vec3 xb = c.free_body ? c.x0n : c.x0;
    f<<"1 Coupled "<<xb(0)<<" "<<xb(1)<<" "<<xb(2)<<" 0 0 0 0 0 0 0 0 0\n";
    f<<"---------------------- POINTS --------------------------------\n";
    f<<"ID Type X Y Z Mass Volume CdA Ca\n";
    f<<"(#) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)\n";
    // MoorDyn points: the free points of the case, then one per line end that is fixed or on the body
    int id = 0;
    std::vector<int> pid(c.points.size());
    for(size_t j=0; j<c.points.size(); ++j)
    {
        const auto &p = c.points[j];
        pid[j] = ++id;
        f<<id<<" Free "<<p.p(0)<<" "<<p.p(1)<<" "<<p.p(2)<<" "<<p.mass<<" "<<p.volume<<" "<<p.CdA<<" "<<p.Ca<<"\n";
    }
    std::vector<std::array<int,2>> le(c.lines.size());
    for(size_t k=0; k<c.lines.size(); ++k)
    {
        const end_ref *r[2] = {&c.lines[k].a, &c.lines[k].b};
        for(int e=0; e<2; ++e)
        {
            if(r[e]->kind==2) {le[k][e] = pid[r[e]->idx]; continue;}
            le[k][e] = ++id;
            const Vec3 &p = r[e]->p;
            f<<id<<" "<<(r[e]->kind==0 ? "Fixed" : "Body1")<<" "<<p(0)<<" "<<p(1)<<" "<<p(2)<<" 0 0 0 0\n";
        }
    }
    f<<"---------------------- LINES ----------------------------------------\n";
    f<<"ID LineType AttachA AttachB UnstrLen NumSegs LineOutputs\n";
    f<<"(#) (name) (#) (#) (m) (-) (-)\n";
    for(size_t k=0; k<c.lines.size(); ++k)
    {
        const auto &l = c.lines[k];
        f<<k+1<<" "<<c.types[l.type].name<<" "<<le[k][0]<<" "<<le[k][1]<<" "<<l.L<<" "<<l.nseg<<" -\n";
    }
    f<<"---------------------- OPTIONS -----------------------------------------\n";
    f<<"0 writeLog\n";
    f<<dtM<<" dtM\n";
    f<<"RK4 tScheme\n";
    f<<c.g<<" g\n";
    f<<c.rho<<" rho\n";
    f<<c.depth<<" WtrDpth\n";
    f<<c.kb<<" kBot\n";
    f<<c.cb<<" cBot\n";
    f<<(c.md_ic_dynamic ? 1 : 0)<<" ICgenDynamic\n";
    f<<"1.0 dtIC\n";
    f<<(c.md_ic_dynamic ? c.md_ic_time : 300.0)<<" TmaxIC\n";
    f<<"5.0 CdScaleIC\n";
    f<<(c.md_ic_dynamic ? "1.0e-9" : "1.0e-6")<<" threshIC\n";
    f<<"0.0 FrictionCoefficient\n";
    if(!c.current.empty())
    {
        f<<"1 Currents\n";
        std::ofstream u(dir + "/current_profile.txt");
        u<<std::setprecision(12)<<"--------------------- MoorDyn steady currents File ----------------------------------\n";
        u<<"current profile of "<<c.name<<"\nz (m), ux (m/s), uy (m/s), uz (m/s)\n";
        for(const auto &r : c.current)
            u<<r[0]<<" "<<r[1]<<" "<<r[2]<<" "<<r[3]<<"\n";
    }
    if(!c.waves.empty())
        f<<"1 WaveKin\n";             // kinematics at the nodes through the API (CORAL's wave field)
    f<<"1 disableOutput\n";
    f<<"1 disableOutTime\n";
    f<<"------------------------- need this line --------------------------------------\n";
}

// ------------------------------------------------------------------ MoorDyn run

struct md_handles
{
    MoorDyn sys = nullptr;
    std::vector<MoorDynLine> lines;
    std::vector<MoorDynPoint> points;
    std::vector<unsigned int> N;
};

static std::vector<double> md_row(const case_spec &c, const md_handles &H, double t, const double *x, const double *f)
{
    std::vector<double> r = {t, x[0], x[2], x[4], f[0], f[1], f[2], f[4], std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2])};
    for(size_t k=0; k<c.lines.size(); ++k)
    {
        double Ta[3], Tb[3];
        MoorDyn_GetLineNodeTen(H.lines[k], 0, Ta);
        MoorDyn_GetLineNodeTen(H.lines[k], H.N[k], Tb);
        r.push_back(std::sqrt(Ta[0]*Ta[0] + Ta[1]*Ta[1] + Ta[2]*Ta[2]));
        r.push_back(std::sqrt(Tb[0]*Tb[0] + Tb[1]*Tb[1] + Tb[2]*Tb[2]));
    }
    for(size_t j=0; j<c.points.size(); ++j)
    {
        double p[3];
        MoorDyn_GetPointPos(H.points[j], p);
        r.push_back(p[0]); r.push_back(p[1]); r.push_back(p[2]);
    }
    return r;
}

static series run_moordyn(const case_spec &c, const std::string &dir, double dtM)
{
    write_moordyn_input(c, dir, dtM);
    md_handles H;
    H.sys = MoorDyn_Create((dir + "/moordyn.txt").c_str());
    if(!H.sys) throw std::runtime_error("MoorDyn_Create failed for " + c.name);
    unsigned int ndof = 0;
    MoorDyn_NCoupledDOF(H.sys, &ndof);
    if(ndof != 6) throw std::runtime_error("MoorDyn: expected 6 coupled DOF, got " + std::to_string(ndof));
    for(size_t k=0; k<c.lines.size(); ++k)
    {
        H.lines.push_back(MoorDyn_GetLine(H.sys, (unsigned int)k + 1));
        unsigned int n;
        MoorDyn_GetLineN(H.lines.back(), &n);
        H.N.push_back(n);
    }
    for(size_t j=0; j<c.points.size(); ++j)
        H.points.push_back(MoorDyn_GetPoint(H.sys, (unsigned int)j + 1));

    double x[6] = {0}, xd[6] = {0}, f[6] = {0};
    // state: prescribed motion, or the free body (x, z, pitch and their rates)
    Vec3 xp, vp, ap;
    double th, thd, thdd;
    double q[3] = {c.x0(0), c.x0(2), 0.0}, qd[3] = {0.0, 0.0, 0.0};
    auto set_x = [&](double t)
    {
        if(c.free_body)
        {
            x[0] = q[0]; x[1] = c.x0(1); x[2] = q[1]; x[3] = 0.0; x[4] = q[2]; x[5] = 0.0;
            xd[0] = qd[0]; xd[1] = 0.0; xd[2] = qd[1]; xd[3] = 0.0; xd[4] = qd[2]; xd[5] = 0.0;
        }
        else
        {
            motion(c, t, xp, vp, ap, th, thd, thdd);
            x[0] = xp(0); x[1] = xp(1); x[2] = xp(2); x[3] = 0.0; x[4] = th; x[5] = 0.0;
            xd[0] = vp(0); xd[1] = vp(1); xd[2] = vp(2); xd[3] = 0.0; xd[4] = thd; xd[5] = 0.0;
        }
    };
    set_x(0.0);
    if(MoorDyn_Init(H.sys, x, xd) != MOORDYN_SUCCESS) throw std::runtime_error("MoorDyn_Init failed for " + c.name);
    // waves: CORAL's linear wave field at MoorDyn's nodes, at the middle of each coupling step
    const auto waves = make_waves(c);
    unsigned int nwp = 0;
    std::vector<double> rwp, uwp, awp;
    if(waves)
    {
        if(MoorDyn_ExternalWaveKinInit(H.sys, &nwp) != MOORDYN_SUCCESS) throw std::runtime_error("MoorDyn: no external wave kinematics");
        rwp.resize(3*nwp); uwp.resize(3*nwp); awp.resize(3*nwp);
    }
    auto set_waves = [&](double tm)
    {
        if(!waves) return;
        MoorDyn_ExternalWaveKinGetCoordinates(H.sys, rwp.data());
        for(unsigned int i=0; i<nwp; ++i)
        {
            Vec3 u, du;
            double pr;
            waves->kinematics(Vec3(rwp[3*i], rwp[3*i+1], rwp[3*i+2]), tm, u, du, pr);
            for(int d=0; d<3; ++d) {uwp[3*i+d] = u(d); awp[3*i+d] = du(d);}
        }
        MoorDyn_ExternalWaveKinSet(H.sys, uwp.data(), awp.data(), tm);
    };
    double t = 0.0, dt0 = 0.0;
    MoorDyn_Step(H.sys, x, xd, f, &t, &dt0);      // dt 0: the forces of the initial state
    series s;
    s.names = channel_names(c);
    s.rows.push_back(md_row(c, H, 0.0, x, f));
    const double h = c.md_h > 0.0 ? c.md_h : c.h;
    const int n = (int)std::lround(c.tend/h);
    const int every = std::max(1, (int)std::lround(c.out_dt/h));
    const double w = 2.0*M_PI/c.F_T;
    for(int i=1; i<=n; ++i)
    {
        double tt = (i - 1)*h, dt = h;
        set_x(tt);
        set_waves(tt + 0.5*h);
        if(MoorDyn_Step(H.sys, x, xd, f, &tt, &dt) != MOORDYN_SUCCESS)
            throw std::runtime_error("MoorDyn_Step failed for " + c.name);
        const double t1 = i*h;
        if(c.free_body)
        {
            // M qdd = F_moor(t1) + F(t1) - K (q - q0) - B qd, semi-implicit Euler; MoorDyn's force at the end of the
            // step, with the ends moved linearly from q, qd over the step
            const double Fx = f[0] + c.F_amp*std::sin(w*t1) - c.K11*(q[0] - c.x0n(0)) - c.B11*qd[0];
            const double Fz = f[2] - c.K33*(q[1] - c.x0n(2)) - c.B33*qd[1];
            const double My = f[4] - c.K55*q[2] - c.B55*qd[2];
            qd[0] += h*Fx/c.mass;
            qd[1] += h*Fz/c.mass;
            qd[2] += h*My/c.Iyy;
            q[0] += h*qd[0];
            q[1] += h*qd[1];
            q[2] += h*qd[2];
            set_x(t1);
        }
        else
            set_x(t1);
        if(i%every==0)
            s.rows.push_back(md_row(c, H, t1, x, f));
    }
    MoorDyn_Close(H.sys);
    return s;
}

int main(int argc, char **argv)
{
    bool ref = false;
    double h_set = 0.0;
    std::vector<std::string> only;
    for(int i=1; i<argc; ++i)
    {
        const std::string a = argv[i];
        if(a=="--ref") ref = true;
        else if(a=="--h" && i+1 < argc) h_set = std::stod(argv[++i]);     // coupling step of all cases
        else only.push_back(a);
    }
    fs::create_directories("out");
    if(ref) fs::create_directories("ref");
    int nfail = 0;
    for(case_spec c : cases())
    {
        if(!only.empty() && std::find(only.begin(), only.end(), c.name)==only.end())
            continue;
        if(h_set > 0.0) c.h = c.md_h = h_set;
        std::cout<<"== "<<c.name<<": "<<c.title<<std::endl;
        double dtl = 0.0;
        const series sc = run_coral(c, &dtl);
        // MoorDyn on the same sub-step (the coupling step a multiple of it)
        const double hm = c.md_h > 0.0 ? c.md_h : c.h;
        const double dtM = hm/std::ceil(hm/dtl - 1.0e-9);
        // MoorDyn's free points start where CORAL's static start put them (its stationary solver stalls from a rough
        // guess; both solve the same equilibrium, which the first rows compare), at most at the seabed: MoorDyn-C
        // lowers the water depth to the lowest point of its input
        case_spec cm = c;
        for(size_t j=0; j<c.points.size(); ++j)
            for(int d=0; d<3; ++d)
                cm.points[j].p(d) = sc.rows[0][sc.col(c.points[j].name + (d==0 ? "_x" : d==1 ? "_y" : "_z"))];
        for(auto &p : cm.points)
            p.p(2) = std::max(p.p(2), -c.depth);
        series sm;
        try {sm = run_moordyn(cm, "out/" + c.name, dtM);}
        catch(const std::exception &e) {std::cout<<"  MoorDyn: "<<e.what()<<std::endl; ++nfail; continue;}
        sc.write("out/" + c.name + "_coral.csv");
        sm.write("out/" + c.name + "_moordyn.csv");
        if(ref && c.unit) sm.write("ref/" + c.name + ".csv", 7);       // the references of the unit test
        std::cout<<"  sub-step "<<dtM<<" s (CORAL "<<dtl<<" s)"<<std::endl;
        double worst = 0.0;
        for(const auto &e : c.peaks ? compare_peaks(sm, sc, c.t_cmp, c.T) : compare(sm, sc, c.t_cmp))
        {
            if(e.name=="x" || e.name=="z" || e.name=="pitch")
                if(!c.free_body) continue;
            std::printf("  %-10s max diff %12.5g  range %12.5g  rel %8.4f  rms %12.5g\n", e.name.c_str(), e.maxdiff, e.range, e.rel, e.rms);
            worst = std::max(worst, e.rel);
        }
        std::printf("  worst %.4f (tolerance of the unit test %.3f)%s\n", worst, c.tol, worst > c.tol ? "  <-- above" : "");
        if(worst > c.tol) ++nfail;
        // published or measured references
        if(!c.ext_ref.empty())
        {
            const auto ref = read_ext(c.ext_ref);
            for(int k=0; k<2; ++k)
            {
                double md, rms, tm;
                compare_ext(k==0 ? sc : sm, ref, c.ext_channel, c.ext_scale, c.ext_t0, md, rms, tm);
                std::printf("  %-8s against %s: max diff %.3f (at %.1f s), rms %.3f, tolerance %.3f%s\n", k==0 ? "CORAL" : "MoorDyn",
                            c.ext_ref.c_str(), md, tm, rms, c.ext_tol, (k==0 && md > c.ext_tol) ? "  <-- above" : "");
            }
        }
        if(c.peaks)
        {
            const double pc = c.ext_scale*mean_peak(sc, c.ext_channel, c.t_cmp, c.T), pm = c.ext_scale*mean_peak(sm, c.ext_channel, c.t_cmp, c.T);
            if(std::isfinite(c.ext_peak))
                std::printf("  mean peak %s: published %.4g, CORAL %.4g (%+.1f %%), MoorDyn %.4g (%+.1f %%)\n", c.ext_channel.c_str(), c.ext_peak,
                            pc, 100.0*(pc/c.ext_peak - 1.0), pm, 100.0*(pm/c.ext_peak - 1.0));
            else
                std::printf("  mean peak %s: CORAL %.4g, MoorDyn %.4g\n", c.ext_channel.c_str(), pc, pm);
        }
        if(std::isfinite(c.ext_static))
        {
            const int ch = sc.col(c.ext_channel);
            const double vc = c.ext_scale*sc.rows[0][ch], vm = c.ext_scale*sm.rows[0][sm.col(c.ext_channel)];
            std::printf("  static %s: measured %.4g, CORAL %.4g (%+.2f %%), MoorDyn %.4g (%+.2f %%)\n", c.ext_channel.c_str(), c.ext_static,
                        vc, 100.0*(vc/c.ext_static - 1.0), vm, 100.0*(vm/c.ext_static - 1.0));
        }
    }
    std::cout<<(nfail ? "cases above tolerance: " + std::to_string(nfail) : std::string("all within tolerance"))<<std::endl;
    return nfail ? 1 : 0;
}
