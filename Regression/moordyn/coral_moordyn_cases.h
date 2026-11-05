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
// REEF3D::CORAL against MoorDyn v2: the cases and the CORAL side of the comparison. The same case runs in
// CORAL (run_coral, here) and in MoorDyn (coral_moordyn_compare.cpp, which also writes the reference time
// series ref/<case>.csv); the unit test Regression/unit/coral_moordyn_test.cpp runs the CORAL side against
// the references without MoorDyn.
//
// A case: line types, free points, lines between fixed points (anchors, global position), the body (ends at
// positions relative to its origin) and points; the body either moves as prescribed (a smooth ramp into a
// harmonic motion, translation and pitch) or is a free planar rigid body (surge, heave, pitch about its
// origin = CoG) with linear stiffness and damping and a harmonic force: then CORAL integrates it with its
// generalized-alpha scheme and the MoorDyn driver with its own scheme on the same equations of motion.
// Still water at z = 0, flat seabed at -depth, an optional current profile (linear in z, held outside).
//
// Channels (one row per output time): t, body x z pitch, the mooring force on the body Fx Fy Fz and the
// moment My about the body origin (with the inertia of the end nodes, as a coupled MoorDyn body), then per
// line the tension of the end segments Ta Tb (with the axial damping), then per point x y z. F: the magnitude of the
// force on the body (the 'top tension' of published results: the end segment and the end half-segment's loads).

#ifndef CORAL_MOORDYN_CASES_H_
#define CORAL_MOORDYN_CASES_H_

#include"coral/core/coral_model.h"
#include"coral/core/coral_timeint.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_load.h"
#include"coral/structure/coral_linear6.h"
#include"coral/mooring/coral_mooring.h"
#include"coral/hydro/coral_waves.h"
#include<limits>
#include<array>
#include<cmath>
#include<fstream>
#include<memory>
#include<sstream>
#include<stdexcept>
#include<string>
#include<vector>

namespace cmd
{
    using coral::Vec3;

    struct line_type
    {
        std::string name;
        double d = 0.1, mass = 0.0, EA = 0.0;
        double BA = 0.0;                     // [N s]; < 0: damping ratio (MoorDyn)
        double EI = 0.0;
        double Cdn = 1.2, Can = 1.0, Cdt = 0.05, Cat = 0.0;
        std::vector<std::array<double,2>> EA_curve;     // (strain, tension) instead of EA
    };

    // an end: kind 0 fixed (p global), 1 the body (p relative to the body origin), 2 a free point (idx)
    struct end_ref
    {
        int kind = 0;
        Vec3 p = Vec3::Zero();
        int idx = -1;
    };

    struct line_spec
    {
        std::string name;
        int type = 0;
        end_ref a, b;
        double L = 0.0;
        int nseg = 20;
    };

    struct point_spec
    {
        std::string name;
        Vec3 p = Vec3::Zero();
        double mass = 0.0, volume = 0.0, CdA = 0.0, Ca = 0.0;
    };

    struct case_spec
    {
        std::string name, title;
        double depth = 100.0, rho = 1025.0, g = 9.81, kb = 3.0e6, cb = 3.0e5;
        std::vector<line_type> types;
        std::vector<point_spec> points;
        std::vector<line_spec> lines;
        std::vector<std::array<double,4>> current;   // (z, ux, uy, uz), increasing z; empty: none
        Vec3 x0 = Vec3::Zero();                      // body origin at t = 0 (after the offset of a free body)

        // prescribed motion: x0 + ramp(t) (amp sin(w t)), pitch ramp(t) pitch_amp sin(w t)
        Vec3 amp = Vec3::Zero();
        double pitch_amp = 0.0, T = 10.0, ramp = 10.0;
        double t_motion = 0.0;                       // the motion (and its ramp) starts at this time
        int ramp_kind = 0;                           // 0: 0.5 (1 - cos(pi t/ramp)), 1: f^3 (6 f^2 - 15 f + 10), f = t/ramp

        // regular wave components (H, T, direction, phase [rad]) in CORAL's convention, eta = H/2 cos(w t - k (x cos b
        // + y sin b) + phase), finite depth; MoorDyn takes CORAL's kinematics at its nodes (WaveKin 1, the free
        // surface stays at z = 0 there, CORAL's nodes see the elevation)
        std::vector<std::array<double,4>> waves;

        // a published or measured reference (CSV 't, value' per line): value = ext_scale * channel ext_channel at
        // t + ext_t0; ext_tol: the largest difference allowed (in the reference's units). ext_static: a measured
        // static value of the channel at t = 0 (NaN: none), ext_static_tol relative
        std::string ext_ref, ext_channel = "L1_Tb";
        double ext_scale = 1.0e-3, ext_t0 = 0.0, ext_tol = 0.0;
        double ext_static = std::numeric_limits<double>::quiet_NaN(), ext_static_tol = 0.02;
        // snap loads (a line that goes slack and snaps taut): the codes are compared by the mean of the largest value
        // of each channel per period of the motion over t >= t_cmp, not by the time series (the instants of the snaps
        // are sensitive); ext_peak: the mean peak of ext_channel expected (published ratio times the quasi-static
        // value), ext_peak_tol relative
        bool peaks = false;
        double ext_peak = std::numeric_limits<double>::quiet_NaN(), ext_peak_tol = 0.2;
        bool unit = true;                            // part of the unit test (else only coral_moordyn_compare)

        // free planar body (surge, heave, pitch): mass and pitch inertia about the origin, stiffness and damping
        // against x0n (neutral), a harmonic force F_amp sin(2 pi t / F_T) at the origin
        bool free_body = false;
        Vec3 x0n = Vec3::Zero();
        double mass = 0.0, Iyy = 0.0, K11 = 0.0, K33 = 0.0, K55 = 0.0, B11 = 0.0, B33 = 0.0, B55 = 0.0;
        double F_amp = 0.0, F_T = 1.0;

        double tend = 60.0, h = 0.01, out_dt = 0.1;
        double dt_line = 0.0;                        // sub-step (0: CORAL's automatic one, also given to MoorDyn)
        double t_cmp = 0.0;                          // compare from this time on
        double tol = 0.02;                           // unit test: largest difference / range of each channel
        // MoorDyn's initial condition: the stationary solver (default), or the dynamic relaxation with the drag
        // scaled by CdScaleIC (md_ic_time: its duration), for lines where the stationary solver stalls
        bool md_ic_dynamic = false;
        double md_ic_time = 300.0;
        // MoorDyn's coupling step (0: h). MoorDyn moves a coupled body linearly over a coupling step from its
        // start velocity (CORAL: cubic Hermite through both ends of the step), an error of a h^2 / 2 at the
        // fairlead that stiff segments turn into tension (EA/l 1e8 N/m: 15 kN at h 0.01 s); the references take 1 ms
        double md_h = 0.001;
    };

    // ------------------------------------------------------------------ time series

    struct series
    {
        std::vector<std::string> names;
        std::vector<std::vector<double>> rows;

        int col(const std::string &n) const
        {
            for(size_t i=0; i<names.size(); ++i) if(names[i]==n) return (int)i;
            return -1;
        }
        void write(const std::string &file, int precision=10) const
        {
            std::ofstream f(file);
            if(!f) throw std::runtime_error("cannot write " + file);
            for(size_t i=0; i<names.size(); ++i) f<<(i ? "," : "")<<names[i];
            f<<"\n";
            f.precision(precision);
            for(const auto &r : rows)
            {
                for(size_t i=0; i<r.size(); ++i) f<<(i ? "," : "")<<r[i];
                f<<"\n";
            }
        }
        static series read(const std::string &file)
        {
            std::ifstream f(file);
            if(!f) throw std::runtime_error("cannot read " + file);
            series s;
            std::string line;
            std::getline(f, line);
            {
                std::stringstream ss(line);
                std::string c;
                while(std::getline(ss, c, ',')) s.names.push_back(c);
            }
            while(std::getline(f, line))
            {
                if(line.empty()) continue;
                std::stringstream ss(line);
                std::string c;
                std::vector<double> r;
                while(std::getline(ss, c, ',')) r.push_back(std::stod(c));
                s.rows.push_back(r);
            }
            return s;
        }
    };

    inline std::vector<std::string> channel_names(const case_spec &c)
    {
        std::vector<std::string> n = {"t", "x", "z", "pitch", "Fx", "Fy", "Fz", "My", "F"};
        for(const auto &l : c.lines) {n.push_back(l.name + "_Ta"); n.push_back(l.name + "_Tb");}
        for(const auto &p : c.points) {n.push_back(p.name + "_x"); n.push_back(p.name + "_y"); n.push_back(p.name + "_z");}
        return n;
    }

    // comparison of two series on the times of a (b interpolated linearly): per channel (not t) the largest
    // difference over t >= t0 relative to the range (max - min) of a there, floored at 'floor' of the largest |a|
    // (a tension that varies little about its mean: relative to a tenth of it);
    // a channel that is zero in a (|a| < 1e-3, Fy of a planar case): the absolute difference must stay below 1e-3
    struct channel_error {std::string name; double maxdiff = 0.0, range = 0.0, rel = 0.0, rms = 0.0;};
    inline std::vector<channel_error> compare(const series &a, const series &b, double t0, double floor=0.1)
    {
        std::vector<channel_error> out;
        for(size_t c=1; c<a.names.size(); ++c)
        {
            const int cb = b.col(a.names[c]);
            if(cb < 0) continue;
            channel_error e;
            e.name = a.names[c];
            double amin = 1e300, amax = -1e300, amag = 0.0, s2 = 0.0;
            int n = 0;
            size_t j = 0;
            for(const auto &ra : a.rows)
            {
                const double t = ra[0];
                if(t < t0 - 1e-9) continue;
                while(j + 1 < b.rows.size() && b.rows[j+1][0] < t) ++j;
                if(j + 1 >= b.rows.size()) break;
                const double tb0 = b.rows[j][0], tb1 = b.rows[j+1][0];
                const double f = tb1 > tb0 ? std::clamp((t - tb0)/(tb1 - tb0), 0.0, 1.0) : 0.0;
                const double vb = b.rows[j][cb] + f*(b.rows[j+1][cb] - b.rows[j][cb]);
                const double d = std::fabs(ra[c] - vb);
                e.maxdiff = std::max(e.maxdiff, d);
                s2 += d*d;
                ++n;
                amin = std::min(amin, ra[c]);
                amax = std::max(amax, ra[c]);
                amag = std::max(amag, std::fabs(ra[c]));
            }
            e.range = std::max(amax - amin, floor*amag);
            e.rms = n ? std::sqrt(s2/n) : 0.0;
            e.rel = amag < 1.0e-3 ? (e.maxdiff < 1.0e-3 ? 0.0 : 1.0e300) : e.maxdiff/e.range;
            out.push_back(e);
        }
        return out;
    }

    // the mean over the periods T (from t0) of the largest value of a channel
    inline double mean_peak(const series &s, const std::string &ch, double t0, double T)
    {
        const int c = s.col(ch);
        double sum = 0.0, mx = -1e300;
        int n = 0, k = -1;
        for(const auto &r : s.rows)
        {
            if(r[0] < t0 - 1e-9) continue;
            const int kp = (int)std::floor((r[0] - t0)/T + 1e-9);
            if(kp != k) {if(k >= 0) {sum += mx; ++n;} k = kp; mx = -1e300;}
            mx = std::max(mx, r[c]);
        }
        return n ? sum/n : 0.0;            // complete periods only
    }

    // snap loads: per channel (the force magnitude and the tensions; a signed component peaks near zero when the line is
    // slack) the relative difference of the mean peaks
    inline std::vector<channel_error> compare_peaks(const series &a, const series &b, double t0, double T)
    {
        std::vector<channel_error> out;
        for(size_t c=1; c<a.names.size(); ++c)
        {
            const std::string &n = a.names[c];
            const bool tension = n=="F" || (n.size() > 3 && (n.compare(n.size()-3, 3, "_Ta")==0 || n.compare(n.size()-3, 3, "_Tb")==0));
            if(b.col(n) < 0 || !tension) continue;
            channel_error e;
            e.name = n;
            const double pa = mean_peak(a, n, t0, T), pb = mean_peak(b, n, t0, T);
            e.maxdiff = std::fabs(pa - pb);
            e.range = std::fabs(pa);
            e.rel = e.range > 1e-12 ? e.maxdiff/e.range : 0.0;
            out.push_back(e);
        }
        return out;
    }

    // ------------------------------------------------------------------ prescribed motion

    inline void ramp(double t, double Tr, double &r, double &rd, double &rdd, int kind=0)
    {
        if(Tr <= 0.0 || t >= Tr) {r = 1.0; rd = rdd = 0.0; return;}
        if(t <= 0.0) {r = rd = rdd = 0.0; return;}
        if(kind==1)
        {
            const double f = t/Tr;
            r = f*f*f*(6.0*f*f - 15.0*f + 10.0);
            rd = 30.0*f*f*(f*f - 2.0*f + 1.0)/Tr;
            rdd = 60.0*f*(2.0*f*f - 3.0*f + 1.0)/(Tr*Tr);
            return;
        }
        const double a = M_PI/Tr;
        r = 0.5*(1.0 - std::cos(a*t));
        rd = 0.5*a*std::sin(a*t);
        rdd = 0.5*a*a*std::cos(a*t);
    }

    // position, velocity, acceleration of the body origin and pitch, its rate and acceleration
    inline void motion(const case_spec &c, double t, Vec3 &x, Vec3 &v, Vec3 &a, double &th, double &thd, double &thdd)
    {
        double r, rd, rdd;
        t -= c.t_motion;
        ramp(t, c.ramp, r, rd, rdd, c.ramp_kind);
        const double w = 2.0*M_PI/c.T, s = std::sin(w*t), co = std::cos(w*t);
        const double f = r*s, fd = rd*s + r*w*co, fdd = rdd*s + 2.0*rd*w*co - r*w*w*s;
        x = c.x0 + f*c.amp;
        v = fd*c.amp;
        a = fdd*c.amp;
        th = f*c.pitch_amp;
        thd = fd*c.pitch_amp;
        thdd = fdd*c.pitch_amp;
    }

    inline coral::Mat3 rot_y(double th)
    {
        coral::Mat3 R;
        R << std::cos(th), 0.0, std::sin(th),
             0.0, 1.0, 0.0,
             -std::sin(th), 0.0, std::cos(th);
        return R;
    }

    // ------------------------------------------------------------------ waves

    inline std::shared_ptr<coral::wave_field> make_waves(const case_spec &c)
    {
        if(c.waves.empty())
            return nullptr;
        auto w = std::make_shared<coral::wave_field>(c.depth, c.g);
        for(const auto &q : c.waves)
            w->add_regular(q[0], q[1], q[2], q[3]);
        return w;
    }

    // ------------------------------------------------------------------ external references

    inline std::vector<std::array<double,2>> read_ext(const std::string &file)
    {
        std::ifstream f(file);
        if(!f) throw std::runtime_error("cannot read " + file);
        std::vector<std::array<double,2>> v;
        std::string l;
        while(std::getline(f, l))
        {
            const size_t h = l.find('#');
            if(h != std::string::npos) l = l.substr(0, h);
            for(char &ch : l) if(ch==',') ch = ' ';
            std::istringstream ls(l);
            double t, x;
            if(ls >> t >> x) v.push_back({t, x});
        }
        return v;
    }

    // the largest and the rms difference of channel ch of s (scaled, at t + t0, linear between rows) to the reference
    inline void compare_ext(const series &s, const std::vector<std::array<double,2>> &ref, const std::string &ch, double scale,
                            double t0, double &maxdiff, double &rms, double &tmax)
    {
        const int c = s.col(ch);
        if(c < 0) throw std::runtime_error("no channel " + ch);
        maxdiff = rms = tmax = 0.0;
        int n = 0;
        size_t j = 0;
        for(const auto &r : ref)
        {
            const double t = r[0] + t0;
            while(j + 1 < s.rows.size() && s.rows[j+1][0] < t) ++j;
            if(j + 1 >= s.rows.size()) break;
            const double ta = s.rows[j][0], tb = s.rows[j+1][0];
            const double f = tb > ta ? std::clamp((t - ta)/(tb - ta), 0.0, 1.0) : 0.0;
            const double v = scale*(s.rows[j][c] + f*(s.rows[j+1][c] - s.rows[j][c]));
            const double d = std::fabs(v - r[1]);
            if(d > maxdiff) {maxdiff = d; tmax = r[0];}
            rms += d*d;
            ++n;
        }
        rms = n ? std::sqrt(rms/n) : 0.0;
    }

    // ------------------------------------------------------------------ CORAL

    inline coral::mooring_props props(const case_spec &c, const line_type &t, double L)
    {
        coral::mooring_props p;
        p.L = L;
        p.d = t.d;
        p.mass = t.mass;
        p.EA = t.EA;
        p.BA = t.BA;
        p.EI = t.EI;
        p.Cdn = t.Cdn; p.Can = t.Can; p.Cdt = t.Cdt; p.Cat = t.Cat;
        p.kb = c.kb; p.cb = c.cb;
        p.rho = c.rho;
        p.BA_slack = true;                     // MoorDyn's axial damping law
        for(const auto &e : t.EA_curve) {p.EA_strain.push_back(e[0]); p.EA_tension.push_back(e[1]);}
        return p;
    }

    // harmonic force on the body origin
    class harmonic_force : public coral::load
    {
    public:
        harmonic_force(int b, const Vec3 &F, double T) : b_(b), F_(F), w_(2.0*M_PI/T) {}
        std::string name() const override {return "harmonic";}
        coral::coupling_kind coupling() const override {return coral::coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
        void eval(const coral::model &m, double t, coral::VecX &f) const override
        {
            m.add_force(f, b_, m.body(b_).x, F_*std::sin(w_*t));
        }
    private:
        int b_;
        Vec3 F_;
        double w_;
    };

    struct coral_run
    {
        coral::model m;
        coral::mooring_system *sys = nullptr;
        std::unique_ptr<coral::mooring_system> own;     // prescribed motion: not a load of the model
        std::unique_ptr<coral::timeint> ti;
        std::vector<int> body_end;                       // per line: the end on the body (0, 1) or -1
    };

    inline void setup_coral(const case_spec &c, coral_run &R)
    {
        coral::model &m = R.m;
        m.gravity = Vec3(0.0, 0.0, -c.g);
        coral::rigid_body b;
        b.mass = c.free_body ? c.mass : 1.0;
        b.inertia = coral::Mat3::Identity();
        if(c.free_body) b.inertia *= c.Iyy;
        b.x = c.free_body ? c.x0n : c.x0;
        const int ib = m.add_body(b);
        auto sys = std::make_unique<coral::mooring_system>(m, make_waves(c), Vec3::Zero(), "moor");
        if(!c.current.empty())
        {
            std::vector<double> z;
            std::vector<Vec3> u;
            for(const auto &r : c.current) {z.push_back(r[0]); u.push_back(Vec3(r[1], r[2], r[3]));}
            sys->set_current_profile(z, u);
        }
        coral::mooring_system &S = *sys;
        // the lines attach at the neutral position (free body) or the start (prescribed), then the body moves
        const Vec3 xb = m.bodies[ib].x;
        for(const auto &p : c.points)
        {
            coral::mooring_point_props pp;
            pp.mass = p.mass;
            pp.volume = p.volume;
            pp.CdA = p.CdA;
            pp.Ca = p.Ca;
            pp.radius = 0.0;                   // MoorDyn's points do not touch the seabed
            S.add_point(p.name, p.p, pp, -c.depth);
        }
        for(const auto &l : c.lines)
        {
            coral::mooring_end e[2];
            const end_ref *r[2] = {&l.a, &l.b};
            int be = -1;
            for(int k=0; k<2; ++k)
            {
                if(r[k]->kind==0) {e[k].body = -1; e[k].p = r[k]->p;}
                else if(r[k]->kind==1) {e[k].body = ib; e[k].p = xb + r[k]->p; be = k;}
                else e[k].point = r[k]->idx;
            }
            S.add_line(m, l.name, e[0], e[1], props(c, c.types[l.type], l.L), l.nseg, -c.depth);
            R.body_end.push_back(be);
        }
        if(c.dt_line > 0.0)
            S.set_dt_line(c.dt_line);
        if(c.free_body)
        {
            m.add_constraint(std::make_unique<coral::joint>(m, coral::joint_type::planar, -1, ib, xb, Vec3::UnitY()));
            coral::Mat6 K = coral::Mat6::Zero(), C = coral::Mat6::Zero();
            K(0,0) = c.K11; K(2,2) = c.K33; K(4,4) = c.K55;
            C(0,0) = c.B11; C(2,2) = c.B33; C(4,4) = c.B55;
            m.add_load(std::make_unique<coral::linear6>(m, ib, xb, K, C));
            if(c.F_amp != 0.0)
                m.add_load(std::make_unique<harmonic_force>(ib, Vec3(c.F_amp, 0.0, 0.0), c.F_T));
            R.sys = static_cast<coral::mooring_system*>(&m.add_load(std::move(sys)));
            // the offset start: the lines in their static equilibrium there (initialize), the body at rest
            m.bodies[ib].x = c.x0;
            coral::timeint_settings ts;
            ts.rho_inf = 0.9;
            R.ti = std::make_unique<coral::timeint>(ts);
            R.ti->initialize(m, 0.0);
        }
        else
        {
            R.own = std::move(sys);
            R.sys = R.own.get();
            Vec3 x, v, a;
            double th, thd, thdd;
            motion(c, 0.0, x, v, a, th, thd, thdd);
            m.bodies[ib].x = x;
            m.bodies[ib].R = rot_y(th);
            R.sys->initialize(m, 0.0);
        }
    }

    inline void set_prescribed(const case_spec &c, coral::model &m, double t)
    {
        Vec3 x, v, a;
        double th, thd, thdd;
        motion(c, t, x, v, a, th, thd, thdd);
        coral::rigid_body &B = m.bodies[0];
        B.x = x; B.v = v; B.vd = a;
        B.R = rot_y(th);
        B.w = Vec3(0.0, thd, 0.0);            // body frame = global for a rotation about y
        B.wd = Vec3(0.0, thdd, 0.0);
    }

    inline std::vector<double> coral_row(const case_spec &c, const coral_run &R, double t)
    {
        const coral::rigid_body &B = R.m.bodies[0];
        std::vector<double> r = {t, B.x(0), B.x(2), std::atan2(B.R(0,2), B.R(0,0))};
        Vec3 F = Vec3::Zero(), M = Vec3::Zero();
        for(size_t k=0; k<c.lines.size(); ++k)
        {
            const int e = R.body_end[k];
            if(e < 0) continue;
            Vec3 Fa, Fb;
            R.sys->line_end_forces((int)k, Fa, Fb);
            const Vec3 Fe = e==0 ? Fa : Fb;
            const auto nodes = R.sys->line_nodes((int)k);
            const Vec3 p = e==0 ? nodes.front() : nodes.back();
            F += Fe;
            M += (p - B.x).cross(Fe);
        }
        r.push_back(F(0)); r.push_back(F(1)); r.push_back(F(2)); r.push_back(M(1)); r.push_back(F.norm());
        for(size_t k=0; k<c.lines.size(); ++k)
        {
            const int N = R.sys->line_nseg((int)k);
            r.push_back(std::fabs(R.sys->line_tension((int)k, 0)));
            r.push_back(std::fabs(R.sys->line_tension((int)k, N-1)));
        }
        for(size_t j=0; j<c.points.size(); ++j)
        {
            const Vec3 p = R.sys->point_position((int)j);
            r.push_back(p(0)); r.push_back(p(1)); r.push_back(p(2));
        }
        return r;
    }

    // the CORAL run of a case: one row per output time
    inline series run_coral(const case_spec &c, double *dt_line=nullptr)
    {
        coral_run R;
        setup_coral(c, R);
        if(dt_line) *dt_line = R.sys->dt_line();
        series s;
        s.names = channel_names(c);
        const int n = (int)std::lround(c.tend/c.h);
        const int every = std::max(1, (int)std::lround(c.out_dt/c.h));
        s.rows.push_back(coral_row(c, R, 0.0));
        for(int i=1; i<=n; ++i)
        {
            const double t = i*c.h;
            if(c.free_body)
                R.ti->step(R.m, t - c.h, c.h);
            else
            {
                set_prescribed(c, R.m, t);
                R.sys->commit(R.m, t, c.h);
            }
            if(i%every==0)
                s.rows.push_back(coral_row(c, R, t));
        }
        return s;
    }

    // ------------------------------------------------------------------ the cases

    // OC3-Hywind chain (Jonkman 2010): 320 m water depth, anchors at 853.87 m, fairleads at 5.2 m radius and 70 m below
    // the surface, 902.2 m of line
    inline line_type oc3_chain()
    {
        line_type t;
        t.name = "oc3";
        t.d = 0.09; t.mass = 77.7066; t.EA = 384.243e6; t.BA = 2.0e6;
        t.Cdn = 1.6; t.Can = 1.0; t.Cdt = 0.1; t.Cat = 0.0;
        return t;
    }

    inline line_type studless_chain()
    {
        line_type t;
        t.name = "chain";
        t.d = 0.252; t.mass = 390.0; t.EA = 1.674e9; t.BA = 8.0e6;
        t.Cdn = 1.37; t.Can = 1.0; t.Cdt = 0.64; t.Cat = 0.0;
        return t;
    }

    inline line_type polyester()
    {
        line_type t;
        t.name = "poly";
        t.d = 0.16; t.mass = 26.0; t.EA = 1.2e8; t.BA = 4.0e5;
        t.Cdn = 1.2; t.Can = 1.0; t.Cdt = 0.008; t.Cat = 0.0;
        return t;
    }

    inline line_type wire()
    {
        line_type t;
        t.name = "wire";
        t.d = 0.08; t.mass = 32.0; t.EA = 5.5e8; t.BA = 4.0e5;
        t.Cdn = 1.2; t.Can = 1.0; t.Cdt = 0.008; t.Cat = 0.0;
        return t;
    }

    inline end_ref fixed(double x, double y, double z) {end_ref e; e.kind = 0; e.p = Vec3(x, y, z); return e;}
    inline end_ref body(double x, double y, double z) {end_ref e; e.kind = 1; e.p = Vec3(x, y, z); return e;}
    inline end_ref point(int j) {end_ref e; e.kind = 2; e.idx = j; return e;}

    inline point_spec pnt(const std::string &n, double x, double y, double z, double m, double V, double CdA, double Ca)
    {
        point_spec p;
        p.name = n; p.p = Vec3(x, y, z); p.mass = m; p.volume = V; p.CdA = CdA; p.Ca = Ca;
        return p;
    }

    inline std::vector<case_spec> cases()
    {
        std::vector<case_spec> v;
        // ---------------- single lines, prescribed motion
        // OC3 line: static, surge, heave
        {
            case_spec c;
            c.name = "oc3_static"; c.title = "OC3 chain at rest (static start)";
            c.depth = 320.0;
            c.types = {oc3_chain()};
            c.lines = {{"L1", 0, fixed(-853.87, 0.0, -320.0), body(-5.2, 0.0, -70.0), 902.2, 20}};
            c.tend = 20.0; c.t_cmp = 0.0; c.tol = 0.002;
            v.push_back(c);
            c.name = "oc3_surge"; c.title = "OC3 chain, fairlead surge 8 m, T 20 s";
            c.amp = Vec3(8.0, 0.0, 0.0); c.T = 20.0; c.ramp = 20.0; c.tend = 100.0; c.t_cmp = 0.0; c.tol = 0.005;
            v.push_back(c);
            c.name = "oc3_heave"; c.title = "OC3 chain, fairlead heave 3 m, T 10 s";
            c.amp = Vec3(0.0, 0.0, 3.0); c.T = 10.0; c.ramp = 10.0; c.tend = 60.0;
            v.push_back(c);
            // MoorDyn's BA/-zeta: the damping ratio 0.8 of the segments' axial mode (the OC3 input of MoorDyn)
            c.name = "oc3_zeta"; c.title = "OC3 chain with BA = -0.8 (damping ratio), fairlead heave 3 m, T 10 s";
            c.types[0].BA = -0.8;
            v.push_back(c);
        }
        // taut polyester: surge with pitch; the same with a stiffening stress-strain curve
        {
            case_spec c;
            c.name = "taut_poly"; c.title = "taut polyester, surge 4 m and pitch 3 deg, T 12 s";
            c.depth = 100.0;
            c.types = {polyester()};
            c.lines = {{"L1", 0, fixed(-250.0, 0.0, -100.0), body(-5.0, 0.0, -10.0), 258.0, 30}};
            c.amp = Vec3(4.0, 0.0, 0.0); c.pitch_amp = 3.0*M_PI/180.0; c.T = 12.0; c.ramp = 12.0; c.tend = 72.0; c.tol = 0.005;
            v.push_back(c);
            c.name = "poly_curve"; c.title = "taut polyester with a stress-strain curve (stiffening), surge 6 m, T 12 s";
            c.types[0].EA_curve = {{0.0, 0.0}, {0.005, 0.35e6}, {0.01, 0.9e6}, {0.02, 2.4e6}, {0.03, 4.3e6}, {0.05, 8.5e6}};
            c.amp = Vec3(6.0, 0.0, 0.0); c.pitch_amp = 0.0;
            v.push_back(c);
        }
        // heavy chain in 50 m, large heave: slack segments and snap loads near the fairlead (CORAL's seabed contact
        // only pushes, MoorDyn's damping also pulls a node rising from the seabed: 1.7 %, 0.17 % with MoorDyn's law)
        {
            case_spec c;
            c.name = "chain_snap"; c.title = "studless chain in 50 m, fairlead heave 4 m, T 6 s";
            c.depth = 50.0;
            c.types = {studless_chain()};
            c.lines = {{"L1", 0, fixed(-400.0, 0.0, -50.0), body(0.0, 0.0, -5.0), 410.0, 41}};
            c.amp = Vec3(0.0, 0.0, 4.0); c.T = 6.0; c.ramp = 6.0; c.tend = 42.0; c.tol = 0.025;
            c.md_ic_dynamic = true;
            v.push_back(c);
        }
        // a chain hanging through the free surface: fairlead 5 m above the water, heave 3 m
        {
            case_spec c;
            c.name = "chain_surface"; c.title = "studless chain through the free surface, fairlead at +5 m, heave 3 m, T 8 s";
            c.depth = 50.0;
            c.types = {studless_chain()};
            c.lines = {{"L1", 0, fixed(-400.0, 0.0, -50.0), body(0.0, 0.0, 5.0), 412.0, 41}};
            c.amp = Vec3(0.0, 0.0, 3.0); c.T = 8.0; c.ramp = 8.0; c.tend = 48.0; c.tol = 0.012;
            c.md_ic_dynamic = true;
            v.push_back(c);
        }
        // OC3 line in a uniform cross current of 1 m/s and in a sheared current along the line; MoorDyn's initial
        // condition leaves the current out (the line drifts into it), CORAL's static start takes it: the steady
        // state, then the motion
        {
            case_spec c;
            c.name = "oc3_current"; c.title = "OC3 chain in a uniform cross current of 1 m/s";
            c.depth = 320.0;
            c.types = {oc3_chain()};
            c.lines = {{"L1", 0, fixed(-853.87, 0.0, -320.0), body(-5.2, 0.0, -70.0), 902.2, 20}};
            c.current = {{-320.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 1.0, 0.0}};
            c.tend = 200.0; c.t_cmp = 180.0; c.tol = 0.002;
            v.push_back(c);
            c.name = "oc3_shear"; c.title = "OC3 chain in a sheared current (0.2 m/s at the seabed, 1.2 m/s at the surface), surge 5 m";
            c.current = {{-320.0, 0.2, 0.0, 0.0}, {-100.0, 0.6, 0.0, 0.0}, {0.0, 1.2, 0.0, 0.0}};
            c.amp = Vec3(5.0, 0.0, 0.0); c.T = 15.0; c.ramp = 15.0; c.t_motion = 180.0; c.tend = 270.0; c.t_cmp = 175.0;
            c.tol = 0.005;
            v.push_back(c);
        }
        // ---------------- line networks, prescribed motion
        // bridle: two wires from the body to a connection point, one chain to the anchor
        {
            case_spec c;
            c.name = "bridle"; c.title = "bridle: two wires to a connection point, chain to the anchor; surge 6 m and pitch 2 deg, T 15 s";
            c.depth = 150.0;
            c.types = {oc3_chain(), wire()};
            c.points = {pnt("P1", -60.0, 0.0, -40.0, 0.0, 0.0, 0.0, 0.0)};
            c.lines = {{"B1", 1, body(-10.0, 8.0, -15.0), point(0), 55.0, 12},
                       {"B2", 1, body(-10.0, -8.0, -15.0), point(0), 55.0, 12},
                       {"M1", 0, point(0), fixed(-560.0, 0.0, -150.0), 520.0, 30}};
            c.amp = Vec3(6.0, 0.0, 0.0); c.pitch_amp = 2.0*M_PI/180.0; c.T = 15.0; c.ramp = 15.0; c.tend = 90.0; c.tol = 0.005;
            v.push_back(c);
        }
        // clump weight between two chains, resting on the seabed (the seabed contact of the end nodes carries it; the
        // anchor chain lifts off and lands: MoorDyn's seabed damping also pulls, 1.7 % here, 0.15 % with its law)
        {
            case_spec c;
            c.name = "clump"; c.title = "clump weight (8 t, 1 m^3) between two chains; heave 2 m and surge 3 m, T 10 s";
            c.depth = 80.0;
            c.types = {oc3_chain()};
            c.points = {pnt("C1", -150.0, 0.0, -60.0, 8000.0, 1.0, 1.5, 0.8)};
            c.lines = {{"A1", 0, fixed(-320.0, 0.0, -80.0), point(0), 190.0, 20},
                       {"F1", 0, point(0), body(0.0, 0.0, -10.0), 160.0, 20}};
            c.amp = Vec3(3.0, 0.0, 2.0); c.T = 10.0; c.ramp = 10.0; c.tend = 70.0; c.tol = 0.025;
            v.push_back(c);
        }
        // a submerged buoy between two chains (MoorDyn's points always take their whole volume: the buoy must stay
        // under water); short axial ringing of the anchor chain where it lifts off the seabed
        {
            case_spec c;
            c.name = "buoy"; c.title = "submerged buoy (6 m^3, CdA 2.5, Ca 1) between two chains; surge 5 m, T 12 s";
            c.depth = 100.0;
            c.types = {oc3_chain()};
            c.points = {pnt("U1", -120.0, 0.0, -45.0, 500.0, 6.0, 2.5, 1.0)};
            c.lines = {{"A1", 0, fixed(-260.0, 0.0, -100.0), point(0), 175.0, 25},
                       {"F1", 0, point(0), body(0.0, 0.0, -20.0), 140.0, 20}};
            c.amp = Vec3(5.0, 0.0, 0.0); c.T = 12.0; c.ramp = 12.0; c.tend = 72.0; c.tol = 0.035;
            v.push_back(c);
        }
        // chain-wire-chain joined at points without mass (CORAL's line sections)
        {
            case_spec c;
            c.name = "chain_wire_chain"; c.title = "chain-wire-chain joined at massless points; surge 6 m, T 20 s";
            c.depth = 200.0;
            c.types = {oc3_chain(), wire()};
            c.points = {pnt("J1", -700.0, 0.0, -195.0, 0.0, 0.0, 0.0, 0.0), pnt("J2", -60.0, 0.0, -60.0, 0.0, 0.0, 0.0, 0.0)};
            c.lines = {{"S1", 0, fixed(-800.0, 0.0, -200.0), point(0), 100.0, 10},
                       {"S2", 1, point(0), point(1), 660.0, 33},
                       {"S3", 0, point(1), body(0.0, 0.0, -15.0), 80.0, 10}};
            c.amp = Vec3(6.0, 0.0, 0.0); c.T = 20.0; c.ramp = 20.0; c.tend = 100.0; c.tol = 0.005;
            v.push_back(c);
        }
        // ---------------- a moored floating body (planar: surge, heave, pitch) with three OC3 lines
        {
            case_spec c;
            c.name = "spar_decay"; c.title = "OC3-like spar on three lines: surge free decay from 10 m";
            c.depth = 320.0;
            c.types = {oc3_chain()};
            const double ang[3] = {M_PI, M_PI/3.0, -M_PI/3.0};
            for(int k=0; k<3; ++k)
            {
                const double ca = std::cos(ang[k]), sa = std::sin(ang[k]);
                c.lines.push_back({"L" + std::to_string(k+1), 0, fixed(853.87*ca, 853.87*sa, -320.0), body(5.2*ca, 5.2*sa, 8.0), 902.2, 20});
            }
            c.free_body = true;
            c.x0n = Vec3(0.0, 0.0, -78.0);
            c.x0 = Vec3(10.0, 0.0, -78.0);
            c.mass = 1.5e7; c.Iyy = 2.5e10;
            c.K33 = 3.32e5; c.K55 = 1.5e9;
            c.B11 = 1.0e5; c.B33 = 1.3e5; c.B55 = 5.0e7;
            c.h = 0.01; c.out_dt = 0.5; c.tend = 300.0; c.tol = 0.005;
            v.push_back(c);
            c.name = "spar_forced"; c.title = "OC3-like spar on three lines: harmonic surge force 0.5 MN, T 30 s";
            c.x0 = c.x0n;
            c.F_amp = 5.0e5; c.F_T = 30.0;
            c.tend = 240.0;
            v.push_back(c);
            c.name = "spar_waves"; c.title = "OC3-like spar on three lines in waves (H 8 m, T 12 s, on the lines), surge force 0.3 MN at T 12 s";
            c.waves = {{8.0, 12.0, 0.0, 0.0}};
            c.F_amp = 3.0e5; c.F_T = 12.0;
            c.tend = 150.0; c.tol = 0.01;
            v.push_back(c);
        }
        // ---------------- waves (CORAL's linear wave field, at MoorDyn's nodes through its API)
        {
            case_spec c;
            c.name = "oc3_waves"; c.title = "OC3 chain in regular waves (H 6 m, T 10 s, along the line), fairlead surge 3 m";
            c.depth = 320.0;
            c.types = {oc3_chain()};
            c.lines = {{"L1", 0, fixed(-853.87, 0.0, -320.0), body(-5.2, 0.0, -70.0), 902.2, 20}};
            c.waves = {{6.0, 10.0, 0.0, 0.0}};
            c.amp = Vec3(3.0, 0.0, 0.0); c.T = 10.0; c.ramp = 10.0; c.tend = 60.0; c.t_cmp = 0.5; c.tol = 0.02;
            v.push_back(c);
        }
        {
            case_spec c;
            c.name = "chain_irregular"; c.title = "studless chain in 50 m in five oblique wave components (H 0.8 - 2 m, T 6 - 12 s), heave 1 m";
            c.depth = 50.0;
            c.types = {studless_chain()};
            c.lines = {{"L1", 0, fixed(-400.0, 0.0, -50.0), body(0.0, 0.0, -5.0), 410.0, 41}};
            c.waves = {{2.0, 12.0, 0.3, 0.0}, {1.6, 10.0, 0.1, 1.3}, {1.4, 8.5, -0.2, 2.9}, {1.0, 7.0, 0.5, 4.4}, {0.8, 6.0, -0.4, 0.7}};
            c.amp = Vec3(0.0, 0.0, 1.0); c.T = 9.0; c.ramp = 9.0; c.tend = 60.0; c.t_cmp = 0.5; c.tol = 0.03;
            c.md_ic_dynamic = true;
            v.push_back(c);
        }
        // ---------------- Lowe and Langley (OMAE 2006): a flexible jumper with bending stiffness, a published time domain
        // solution (digitized for MoorDyn's tests, ext/lowe_langley_2006_case<n>.csv): the fairlead moves 10 m at T 27 s
        // along x, y, z (cases 1-3, quintic ramp over one period), or waves of 5 m amplitude, T 10 s from 0, 45 and 90
        // degrees pass the fixed jumper (cases 4-6, from t = 0); the reference starts after the ramp (27 s, 20 s)
        {
            line_type t;
            t.name = "jumper";
            t.d = 0.396; t.mass = 165.0; t.EA = 500.0e6; t.BA = -0.5; t.EI = 120.8e3;
            t.Cdn = 1.0; t.Can = 1.0; t.Cdt = 0.0; t.Cat = 0.0;
            case_spec c;
            c.depth = 100.0; c.rho = 1000.0; c.g = 9.807;
            c.types = {t};
            c.lines = {{"L1", 0, fixed(-100.0, 0.0, -55.0), body(0.0, 0.0, 0.0), 170.0, 68}};
            c.x0 = Vec3(0.0, 0.0, -5.0);
            c.ext_channel = "F"; c.ext_scale = 1.0e-3; c.ext_tol = 1.0;
            for(int k=1; k<=3; ++k)
            {
                case_spec q = c;
                q.name = "lowe_langley_" + std::to_string(k);
                q.title = std::string("Lowe and Langley 2006 jumper, fairlead motion 10 m, T 27 s along ") + "xyz"[k-1];
                q.amp = Vec3::Zero(); q.amp(k-1) = 10.0;
                q.T = 27.0; q.ramp = 27.0; q.ramp_kind = 1; q.tend = 27.0 + 80.0; q.tol = k==1 ? 0.01 : 0.05;
                q.ext_ref = "ext/lowe_langley_2006_case" + std::to_string(k) + ".csv"; q.ext_t0 = 27.0;
                q.unit = k != 2;
                v.push_back(q);
            }
            const double phase[3] = {230.0, 165.0, 0.0};
            for(int k=4; k<=6; ++k)
            {
                case_spec q = c;
                q.name = "lowe_langley_" + std::to_string(k);
                q.title = "Lowe and Langley 2006 jumper in waves, amplitude 5 m, T 10 s, from " + std::to_string(45*(k-4)) + " degrees";
                // MoorDyn's test: u ~ cos(k x - w t + phi), CORAL: cos(w t - k x + phase), phase = -phi
                q.waves = {{10.0, 10.0, (k-4)*0.25*M_PI, -phase[k-4]*M_PI/180.0}};
                q.tend = 20.0 + 30.0; q.t_cmp = 0.5; q.tol = 0.02;   // Fx, Fy of the jumper vary little: 1 %
                q.ext_ref = "ext/lowe_langley_2006_case" + std::to_string(k) + ".csv"; q.ext_t0 = 20.0;
                q.unit = k==4;
                v.push_back(q);
            }
        }
        // ---------------- Azcona et al. (Ocean Engineering 2017): a 21 m studless chain in the 5 m deep basin of Ecole
        // Centrale de Nantes, fairlead at the water plane moving horizontally in the plane of the line at T 1.58, 3.16,
        // 4.74 s (0.25 m in the paper's table, taken as the stroke: amplitude 0.125 m, see below); anchor 19.364 m
        // (configuration 1) or 19.870 m (2) away. Measured static fairlead tension 8.13 N and 14.48 N; seabed stiffness
        // 20 N/m^2 and damping 0.1 N s/m^2 per length, Cdn 1.4, Cdt 0.67 (equivalent diameter), Can 1, structural damping
        // 0.1 %, fresh water. Of the dynamic tests the text gives ratios only: at T 1.58 s the largest tension is 'between
        // 2 (Configuration 1) and 3 (Configuration 2) times' the quasi-static one and three times larger in configuration 2
        // than in 1. CORAL (and MoorDyn) give 1.76, 2.74 and 3.07 with the amplitude 0.125 m, 3.8, 6.9 and 3.9 with 0.25 m
        {
            line_type t;
            t.name = "chain";
            t.d = 0.0034; t.mass = 0.069; t.EA = 3.4e5; t.BA = -0.001;
            t.Cdn = 1.4; t.Can = 1.0; t.Cdt = 0.67; t.Cat = 0.0;
            const double dist[2] = {19.364, 19.870}, Tst[2] = {8.13, 14.48};
            for(int k=0; k<2; ++k)
            {
                case_spec c;
                c.depth = 5.0; c.rho = 1000.0; c.g = 9.81;
                c.kb = 20.0/t.d; c.cb = 0.1/t.d;
                c.types = {t};
                c.lines = {{"L1", 0, fixed(-dist[k], 0.0, -5.0), body(0.0, 0.0, 0.0), 21.0, 40}};
                c.h = 0.002; c.out_dt = 0.02; c.md_h = 0.001;
                case_spec q = c;
                q.name = "azcona_c" + std::to_string(k+1) + "_static";
                q.title = "Azcona et al. 2017 chain, configuration " + std::to_string(k+1) + ", static";
                q.tend = 2.0; q.tol = 0.01;          // MoorDyn's stationary solver stops 0.5 % short
                q.ext_channel = "L1_Tb"; q.ext_scale = 1.0; q.ext_static = Tst[k]; q.ext_static_tol = 0.03;
                v.push_back(q);
                const double per[3] = {1.58, 3.16, 4.74};
                // the quasi-static peak (fairlead 0.125 m out, mooring_qs) 9.22 N and 18.20 N
                const double qs_peak[2] = {9.22, 18.20}, ratio[2] = {2.0, 3.0};
                for(int j=0; j<3; ++j)
                {
                    q = c;
                    q.name = "azcona_c" + std::to_string(k+1) + "_T" + std::to_string((int)std::lround(100*per[j]));
                    q.title = "Azcona et al. 2017 chain, configuration " + std::to_string(k+1) + ", fairlead stroke 0.25 m at T " + std::to_string(per[j]).substr(0, 4) + " s";
                    q.amp = Vec3(0.125, 0.0, 0.0); q.T = per[j]; q.ramp = 2.0*per[j];
                    q.tend = 25.0*per[j]; q.t_cmp = 15.0*per[j]; q.peaks = true; q.tol = 0.05;
                    q.ext_channel = "L1_Tb"; q.ext_scale = 1.0;
                    if(j==0) {q.ext_peak = ratio[k]*qs_peak[k]; q.ext_peak_tol = 0.2;}
                    q.unit = j==0;
                    v.push_back(q);
                }
            }
        }
        return v;
    }
}

#endif
