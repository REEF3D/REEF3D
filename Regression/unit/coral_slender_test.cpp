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
// REEF3D::CORAL unit test of the slender-body loads for steep waves and ringing (coral/hydro/coral_member_morison,
// model morison | lagrangian | extended | rainey | kf; coral/hydro/coral_waves gradients):
// the velocity gradient and the surface slope of the linear waves against differences of their kinematics (finite
// and infinite depth, two directions); a fixed pile in a deep-water regular wave: the particle acceleration adds
// nothing there (u.grad u = grad |u|^2/2 = 0), the axial-divergence force against its integral, Rainey's point
// force -1/2 rho A Ca u^2 eta_x at z = eta and Kristiansen-Faltinsen's rho A Ca (4/g) u^2 du/dt at z = 0 against
// their formulas, F_psi = 8 F_eta for small amplitudes; the same pile with a sampled fluid (the stencil gradients,
// the waterline points and their Taylor step) against the linear waves' gradients; the keyword and its channels.
// Architect: Hans Bihs
#include"coral/hydro/coral_member_morison.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include"coral/hydro/coral_waves.h"
#include"coral/io/coral_input.h"
#include<cmath>
#include<cstdio>
#include<iostream>
#include<string>
#include<vector>

using namespace coral;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}
static std::string fmt(const char *f, double a, double b=0.0, double c=0.0, double d=0.0)
{
    char s[320];
    std::snprintf(s, sizeof(s), f, a, b, c, d);
    return s;
}
static const double RHO = 1025.0, G = 9.81;

static void test_gradients()
{
    std::cout<<"gradients of the linear waves"<<std::endl;
    double e = 0.0, es = 0.0, ediv = 0.0, ecurl = 0.0;
    for(double depth : {20.0, -1.0})
    {
        wave_field w(depth, G);
        w.add_component(0.8, 2.0*M_PI/7.0, 0.0, 0.3);
        w.add_component(0.5, 2.0*M_PI/5.0, 30.0*M_PI/180.0, 1.1);
        const double h = 1.0e-5;
        for(const Vec3 &p : {Vec3(1.3, -2.0, -3.5), Vec3(-4.0, 0.7, -12.0)})
            for(double t : {0.4, 3.3})
            {
                Mat3 Gr, Gn;
                Vec3 sl, u1, u2, du;
                double pr;
                w.gradients(p, t, Gr, sl);
                for(int j=0; j<3; ++j)
                {
                    w.kinematics(p + h*Vec3::Unit(j), t, u1, du, pr);
                    w.kinematics(p - h*Vec3::Unit(j), t, u2, du, pr);
                    Gn.col(j) = (u1 - u2)/(2.0*h);
                }
                const double sx = (w.elevation(p(0) + h, p(1), t) - w.elevation(p(0) - h, p(1), t))/(2.0*h);
                const double sy = (w.elevation(p(0), p(1) + h, t) - w.elevation(p(0), p(1) - h, t))/(2.0*h);
                e = std::max(e, (Gr - Gn).norm()/Gn.norm());
                es = std::max(es, std::hypot(sl(0) - sx, sl(1) - sy)/std::hypot(sx, sy));
                ediv = std::max(ediv, std::fabs(Gr.trace())/Gr.norm());
                ecurl = std::max(ecurl, (Gr - Gr.transpose()).norm()/Gr.norm());
            }
    }
    check(e < 1e-7 && es < 1e-7 && ediv < 1e-12 && ecurl < 1e-12,
          fmt("velocity gradient and slope = differences of the kinematics (%.1e, %.1e); div u %.0e, curl u %.0e", e, es, ediv, ecurl));
}

static member_morison *find_member(const model &m)
{
    for(const auto &l : m.load_list)
        if(auto *p = dynamic_cast<member_morison*>(l.get()))
            return p;
    return nullptr;
}

// a fixed pile D 6 m from -60 m to +10 m, 1 m elements with 4 strips, in a regular wave
static std::string pile(double H, double T, const std::string &depth, const std::string &model, const std::string &extra="",
                        int top=10)
{
    return "rho 1025\ngravity 0 0 -9.81\nwaves regular H " + std::to_string(H) + " T " + std::to_string(T) + "\n" + depth
           + "section pile tube D 6 t 0.06 E 2.1e11 G 8.08e10 rho 7850\n"
             "member mp p1 0 0 -60 p2 0 0 " + std::to_string(top) + " n " + std::to_string(60 + top)
           + " section pile orient 1 0 0 end1 ground\n"
             "morison wave member mp Cd 0 Ca 1.0 n 4 model " + model + " " + extra + "\n";
}

static void test_deep_water()
{
    std::cout<<"fixed pile in a deep-water regular wave"<<std::endl;
    const double H = 4.0, T = 9.0, A = 0.25*M_PI*36.0, a = 0.5*H, w = 2.0*M_PI/T, k = w*w/G;
    coral_case cm, cl, ce, cr, ck;
    cm.read_string(pile(H, T, "", "morison"));
    cl.read_string(pile(H, T, "", "lagrangian"));
    ce.read_string(pile(H, T, "", "extended"));
    cr.read_string(pile(H, T, "", "rainey"));
    ck.read_string(pile(H, T, "", "kf"));
    member_morison *mm = find_member(cm.m), *ml = find_member(cl.m), *me = find_member(ce.m), *mr = find_member(cr.m), *mk = find_member(ck.m);
    double el = 0.0, ee = 0.0, er = 0.0, ek = 0.0, fa = 0.0, aa = 0.0, ra = 0.0, ka = 0.0;
    for(int i=0; i<72; ++i)
    {
        const double t = 20.0 + T*i/72.0, th = w*t, eta = a*std::cos(th);
        Vec3 Fm, Fl, Fe, Mm, M;
        mm->total(cm.m, t, Fm, Mm);
        ml->total(cl.m, t, Fl, M);
        me->total(ce.m, t, Fe, M);
        el = std::max(el, std::fabs(Fl(0) - Fm(0)));
        fa = std::max(fa, std::fabs(Fm(0)));
        // axial divergence: rho A Ca u w_z, u = a w e^kz cos, w_z = -a w k e^kz sin; z = 0 kinematics above z = 0
        const double I = eta < 0.0 ? (std::exp(2*k*eta) - std::exp(-120*k))/(2*k) : (1.0 - std::exp(-120*k))/(2*k) + eta;
        const double Fax = -RHO*A*a*a*w*w*k*std::cos(th)*std::sin(th)*I;
        ee = std::max(ee, std::fabs(Fe(0) - Fl(0) - Fax));
        aa = std::max(aa, std::fabs(Fax));
        // Rainey: -1/2 rho A u(eta)^2 eta_x;  KF: rho A (4/g) u(0)^2 du/dt(0)
        const double ue = a*w*std::exp(k*std::min(eta, 0.0))*std::cos(th);
        const double Fr = -0.5*RHO*A*ue*ue*a*k*std::sin(th);
        const double Fk = RHO*A*4.0/G*std::pow(a*w*std::cos(th), 2)*(-a*w*w*std::sin(th));
        er = std::max(er, std::fabs(mr->point_force(cr.m, t)(0) - Fr));
        ek = std::max(ek, std::fabs(mk->point_force(ck.m, t)(0) - Fk));
        ra = std::max(ra, std::fabs(Fr));
        ka = std::max(ka, std::fabs(Fk));
    }
    check(el < 1e-9*fa, fmt("lagrangian = morison in deep water (u.grad u = 0): max difference %.1e of %.3e N", el, fa));
    check(ee < 1e-4*aa, fmt("extended: the axial-divergence force against its integral, max error %.1e of %.3e N", ee/aa, aa));
    check(er < 1e-9*ra && ek < 1e-9*ka, fmt("rainey F_eta (amplitude %.0f N) and kf F_psi (%.0f N) = their formulas (errors %.0e, %.0e)", ra, ka, er/ra, ek/ka));
    // small amplitude: F_psi = 8 F_eta
    coral_case sr, sk;
    sr.read_string(pile(0.02, T, "", "rainey"));
    sk.read_string(pile(0.02, T, "", "kf"));
    double mr2 = 0.0, mk2 = 0.0;
    for(int i=0; i<72; ++i)
    {
        const double t = 20.0 + T*i/72.0;
        mr2 = std::max(mr2, std::fabs(find_member(sr.m)->point_force(sr.m, t)(0)));
        mk2 = std::max(mk2, std::fabs(find_member(sk.m)->point_force(sk.m, t)(0)));
    }
    check(std::fabs(mk2/mr2 - 8.0) < 1e-3, fmt("small amplitude (H 2 cm): F_psi/F_eta = %.4f (8 in the limit)", mk2/mr2));
}

// one fluid step of an emulated resolved fluid (the linear waves): above the free surface the surface velocity
static void emulate(const model &m, sampled_kinematics &sk, const wave_field &w, double t)
{
    std::vector<Vec3> pts;
    sk.collect(m, pts);
    std::vector<Vec3> u(pts.size());
    std::vector<double> eta(pts.size());
    std::vector<int> in(pts.size(), 1);
    for(size_t q=0; q<pts.size(); ++q)
    {
        eta[q] = w.elevation(pts[q](0), pts[q](1), t);
        Vec3 p = pts[q], du;
        p(2) = std::min(p(2), eta[q]);
        double pr;
        w.kinematics(p, t, u[q], du, pr);
    }
    sk.add_sample(t, u, eta, in);
}

static void test_sampled()
{
    std::cout<<"sampled fluid: stencil gradients and waterline points"<<std::endl;
    // the emulated fluid takes the linear waves' kinematics of z = 0 above the still water level, where their
    // vertical derivatives differ from the linear waves' gradients: the distributed part on a submerged member
    // (top at -2 m), the point forces with small waves
    const std::string depth = "waves depth 30\n";
    for(const char *mo : {"extended", "rainey", "kf"})
    {
        const bool pt = std::string(mo) != "extended";
        const double H = pt ? 0.1 : 2.0, T = 8.0;
        const int top = pt ? 10 : -2;
        coral_case ca, cb, c0, c1;
        ca.read_string(pile(H, T, depth, mo, "", top));
        cb.read_string(pile(H, T, depth, mo, pt ? "step 0.05" : "step 0.5", top) + "fluid_kinematics order 2\n");
        c0.read_string(pile(H, T, depth, "morison", "", top));
        c1.read_string(pile(H, T, depth, "morison", "", top) + "fluid_kinematics order 2\n");
        member_morison *ma = find_member(ca.m), *mb = find_member(cb.m), *m0 = find_member(c0.m), *m1 = find_member(c1.m);
        std::vector<Vec3> pts;
        cb.fluid_kin->collect(cb.m, pts);
        const int ns = 4*(60 + top), npts = ns*7 + (pt ? 7 : 0);
        const double dt = 0.02;
        double ed = 0.0, nd = 0.0, ep = 0.0, np = 0.0;
        for(int i=0; i<=(int)std::lround(24.0/dt); ++i)
        {
            const double t = i*dt;
            emulate(cb.m, *cb.fluid_kin, *ca.waves, t);
            emulate(c1.m, *c1.fluid_kin, *ca.waves, t);
            if(t < 16.0)
                continue;
            Vec3 Fa, Fb, F0, F1, M;
            ma->total(ca.m, t, Fa, M);
            mb->total(cb.m, t, Fb, M);
            m0->total(c0.m, t, F0, M);
            m1->total(c1.m, t, F1, M);
            if(!pt)
            {
                // the higher-order part of the distributed force: extended - Morison, each with its own kinematics
                ed = std::max(ed, std::fabs((Fb(0) - F1(0)) - (Fa(0) - F0(0))));
                nd = std::max(nd, std::fabs(Fa(0) - F0(0)));
            }
            else
            {
                ep = std::max(ep, std::fabs(mb->point_force(cb.m, t)(0) - ma->point_force(ca.m, t)(0)));
                np = std::max(np, std::fabs(ma->point_force(ca.m, t)(0)));
            }
        }
        if(!pt)
            check((int)pts.size()==npts && ed < 0.01*nd,
                  fmt("extended, submerged member: %.0f points (strips and six neighbours each); the higher-order force %.3e N, error %.1e of it",
                      (double)pts.size(), nd, ed/nd));
        else
            check((int)pts.size()==npts && ep < 0.02*np,
                  std::string(mo) + fmt(": %.0f points with the waterline's seven; point force %.3e N (H 0.1 m), error %.1e of it",
                                        (double)pts.size(), np, ep/np));
    }
}

static void test_keyword()
{
    std::cout<<"keyword"<<std::endl;
    coral_case c;
    c.read_string(pile(4.0, 9.0, "waves depth 40\n", "rainey"));
    std::vector<std::string> n;
    std::vector<double> v;
    find_member(c.m)->channel_names(n);
    find_member(c.m)->channel_values(c.m, 3.0, v);
    check(n.size()==9 && n[6]=="wave_Fpx" && v.size()==9 && std::fabs(v[6]) > 0.0, "model rainey: channels _Fx .. _Mz and the point force _Fpx _Fpy _Fpz");
    bool err = false;
    try {coral_case d; d.read_string(pile(4.0, 9.0, "", "fnv"));}
    catch(const std::exception &e) {err = std::string(e.what()).find("model must be") != std::string::npos;}
    check(err, "an unknown model is an error");
}

int main()
{
    test_gradients();
    test_deep_water();
    test_sampled();
    test_keyword();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
