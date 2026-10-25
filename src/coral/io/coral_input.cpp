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

#include"coral/io/coral_input.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_gravity.h"
#include"coral/core/coral_force.h"
#include"coral/core/coral_rotation.h"
#include"coral/structure/coral_spring.h"
#include"coral/structure/coral_linear6.h"
#include"coral/structure/coral_frame.h"
#include"coral/structure/coral_superelement.h"
#include"coral/soil/coral_soil.h"
#include"coral/hydro/coral_hydro_db.h"
#include"coral/hydro/coral_hydro_load.h"
#include"coral/hydro/coral_morison.h"
#include"coral/hydro/coral_strip_hydro.h"
#include"coral/hydro/coral_member_morison.h"
#include"coral/core/coral_points.h"
#include"coral/pto/coral_pto.h"
#include"coral/pto/coral_joint_drive.h"
#include"coral/pto/coral_joint_force.h"
#include"coral/hydro/coral_fluid_load.h"
#include"coral/hydro/coral_body_drag.h"
#include"coral/hydro/coral_sampled_kinematics.h"
#include"coral/hydro/coral_member_fluid_load.h"
#include"coral/mooring/coral_mooring.h"
#include"coral/mooring/coral_seabed.h"
#include"coral/wind/coral_wind.h"
#include"coral/wind/coral_wind_grid.h"
#include"coral/wind/coral_wind_probe.h"
#include"coral/aero/coral_rotor.h"
#include"coral/aero/coral_rotor_tune.h"
#include"coral/control/coral_turbine_controller.h"
#include<cmath>
#include<fstream>
#include<functional>
#include<limits>
#include<set>
#include<sstream>
#include<stdexcept>

namespace coral
{

namespace
{
    // one input line: keyword, then words and key-value pairs
    struct item
    {
        int line = 0;
        std::vector<std::string> t;
        int first = 1;      // first word searched for keys: 2 after the name of a named item (a name may be a key)

        [[noreturn]] void fail(const std::string &msg) const
        {
            throw std::runtime_error("CORAL input, line " + std::to_string(line) + " (" + t[0] + "): " + msg);
        }
        int find(const std::string &key, int from=-1) const
        {
            for(int i=(from < 0 ? first : from); i<(int)t.size(); ++i)
                if(t[i]==key)
                    return i;
            return -1;
        }
        bool has(const std::string &key) const {return find(key) >= 0;}
        double num(int i) const
        {
            if(i >= (int)t.size()) fail("value missing at the end of the line");
            try {size_t n; const double v = std::stod(t[i], &n); if(n != t[i].size()) throw 0; return v;}
            catch(...) {fail("'" + t[i] + "' is not a number");}
        }
        double get(const std::string &key, double def) const
        {
            const int i = find(key);
            return i < 0 ? def : num(i+1);
        }
        double need(const std::string &key) const
        {
            const int i = find(key);
            if(i < 0) fail("'" + key + "' is required");
            return num(i+1);
        }
        Vec3 vec(int i) const {return Vec3(num(i), num(i+1), num(i+2));}
        Vec3 vec(const std::string &key, const Vec3 &def) const
        {
            const int i = find(key);
            return i < 0 ? def : vec(i+1);
        }
        Vec3 need_vec(const std::string &key) const
        {
            const int i = find(key);
            if(i < 0) fail("'" + key + "' is required");
            return vec(i+1);
        }
        std::string word(int i) const
        {
            if(i >= (int)t.size()) fail("word missing at the end of the line");
            return t[i];
        }
        std::string word(const std::string &key) const
        {
            const int i = find(key);
            if(i < 0) fail("'" + key + "' is required");
            return word(i+1);
        }
    };
}

void coral_case::read(const std::string &file)
{
    std::ifstream in(file);
    if(!in)
        throw std::runtime_error("CORAL input: cannot open " + file);
    std::stringstream ss;
    ss << in.rdbuf();
    const size_t slash = file.find_last_of('/');
    read_string(ss.str(), slash==std::string::npos ? "." : file.substr(0, slash));
}

void coral_case::read_string(const std::string &text, const std::string &dir)
{
    std::vector<std::pair<int, std::vector<std::string>>> lines;
    std::istringstream in(text);
    std::string raw, acc;
    int ln = 0, start = 0;
    while(std::getline(in, raw))
    {
        ++ln;
        const size_t c = raw.find('#');
        if(c != std::string::npos) raw = raw.substr(0, c);
        while(!raw.empty() && (raw.back()=='\r' || raw.back()==' ' || raw.back()=='\t')) raw.pop_back();
        if(acc.empty()) start = ln;
        bool cont = false;
        if(!raw.empty() && raw.back()=='\\') {raw.pop_back(); cont = true;}
        acc += " " + raw;
        if(cont) continue;
        std::istringstream ws(acc);
        std::vector<std::string> t;
        std::string w;
        while(ws >> w) t.push_back(w);
        if(!t.empty()) lines.push_back({start, t});
        acc.clear();
    }
    parse(lines, dir);
}

void coral_case::parse(const std::vector<std::pair<int, std::vector<std::string>>> &lines, const std::string &dir)
{
    std::vector<item> items;
    // items whose second word is a name (of a body, joint, load, ...): keys are searched after it
    static const std::set<std::string> named = {"body", "joint", "pto", "drive", "jforce", "spring", "linear6", "force", "drag",
        "mooring", "morison", "section", "member", "superelement", "pile", "point_mass", "fluid_body", "fluid_member", "hydro",
        "point", "line_type", "wind_probe", "rotor", "controller", "strip_member"};
    for(const auto &l : lines)
    {
        items.push_back(item{l.first, l.second});
        if(!l.second.empty() && named.count(l.second[0]))
            items.back().first = 2;
    }
    auto path = [&](const std::string &f) {return (f.empty() || f[0]=='/') ? f : dir + "/" + f;};
    auto body_of = [&](const item &it, const std::string &name) -> int
    {
        if(name=="ground") return -1;
        const int b = m.find_body(name);
        if(b < 0) it.fail("unknown body '" + name + "'");
        return b;
    };
    const double deg = M_PI/180.0;

    // pass 1: settings, waves, bodies
    std::vector<std::pair<int, Vec3>> offsets;    // initial displacement of bodies, applied after the set-up
    std::vector<std::pair<mooring_system*, double>> settles;
    std::map<std::string, beam_section> sections;
    bool gravity_on = true;
    int fluid_kin_order = -1;            // fluid_kinematics: order of the time polynomial, -1: not used
    double wdepth = -1.0, wramp = 0.0;
    std::vector<item> wave_items, wind_items;
    for(const item &it : items)
    {
        const std::string &k = it.t[0];
        if(k=="wind")
        {
            if(it.t.size() < 2) it.fail("wind type missing");
            wind_items.push_back(it);
        }
        else if(k=="rho_air") rho_air = it.num(1);
        else if(k=="gravity")
        {
            if(it.t.size()>1 && it.t[1]=="off") {gravity_on = false; m.gravity = Vec3::Zero();}
            else m.gravity = it.vec(1);
        }
        else if(k=="rho") rho = it.num(1);
        else if(k=="time")
        {
            run.dt = it.need("dt");
            run.t_end = it.need("end");
            run.t_start = it.get("start", 0.0);
            run.ti.rho_inf = it.get("rho_inf", run.ti.rho_inf);
            run.ti.tol = it.get("tol", run.ti.tol);
            run.ti.full_newton = it.has("full_newton");
        }
        else if(k=="output")
        {
            if(it.t.size() > 1 && it.t[1]=="vtp")
                run.vtp_dt = it.num(2);
            else
            {
                run.out_dt = it.get("dt", 0.0);
                if(it.has("file")) run.out_file = it.word("file");
            }
        }
        else if(k=="static")
            run.static_ini = true;
        else if(k=="modes")
        {
            run.modes = (int)it.num(1);
            run.modes_damped = it.has("damped");
        }
        else if(k=="section")
        {
            const std::string name = it.word(1);
            beam_section s;
            if(it.t.size() > 2 && it.t[2]=="tube")
                s = beam_section::tube(it.need("D"), it.need("t"), it.need("E"), it.need("G"), it.need("rho"));
            else
            {
                s.EA = it.need("EA"); s.EIy = it.need("EIy"); s.EIz = it.need("EIz"); s.GJ = it.need("GJ"); s.m = it.need("m");
                s.GAy = it.get("GAy", 0.0); s.GAz = it.get("GAz", 0.0);
                s.rIx = it.get("rIx", 0.0); s.rIy = it.get("rIy", 0.0); s.rIz = it.get("rIz", 0.0); s.D = it.get("D", 0.0);
            }
            sections[name] = s;
        }
        else if(k=="fluid_kinematics")
        {
            fluid_kin_order = (int)it.get("order", 2.0);
            if(fluid_kin_order < 0 || fluid_kin_order > 2) it.fail("order must be 0, 1 or 2");
        }
        else if(k=="waves")
        {
            if(it.t.size() < 2) it.fail("wave type missing");
            if(it.t[1]=="depth") wdepth = it.num(2);
            else if(it.t[1]=="ramp") wramp = it.num(2);
            else wave_items.push_back(it);
        }
        else if(k=="body")
        {
            rigid_body b;
            b.name = it.word(1);
            if(m.find_body(b.name) >= 0) it.fail("body '" + b.name + "' defined twice");
            b.mass = it.need("mass");
            const int ii = it.find("inertia");
            if(ii >= 0)
            {
                b.inertia = it.vec(ii+1).asDiagonal();
                if(ii+6 < (int)it.t.size())
                {
                    bool more = true;
                    try {std::stod(it.t[ii+4]);} catch(...) {more = false;}
                    if(more)
                    {
                        const Vec3 o = it.vec(ii+4);
                        b.inertia(0,1) = b.inertia(1,0) = o(0);
                        b.inertia(0,2) = b.inertia(2,0) = o(1);
                        b.inertia(1,2) = b.inertia(2,1) = o(2);
                    }
                }
            }
            b.x = it.vec("cog", Vec3::Zero());
            b.R = matrix_from_euler(it.vec("rpy", Vec3::Zero())*deg);
            b.v = it.vec("velocity", Vec3::Zero());
            b.w = it.vec("omega", Vec3::Zero());
            const int bi = m.add_body(b);
            if(it.has("offset"))
                offsets.push_back({bi, it.need_vec("offset")});
        }
    }
    if(!wave_items.empty())
    {
        waves = std::make_shared<wave_field>(wdepth, std::fabs(m.gravity(2)) > 0.0 ? std::fabs(m.gravity(2)) : 9.81);
        waves->set_ramp(wramp);
        for(const item &it : wave_items)
        {
            const std::string &type = it.t[1];
            const double dir = it.get("direction", 0.0)*deg;
            if(type=="regular")
                waves->add_regular(it.need("H"), it.need("T"), dir, it.get("phase", 0.0)*deg);
            else if(type=="component")
                waves->add_component(it.need("a"), it.need("w"), dir, it.get("phase", 0.0)*deg);
            else if(type=="jonswap" || type=="pm")
            {
                const double Tp = it.need("Tp");
                const double wp = 2.0*M_PI/Tp;
                waves->add_spectrum(type=="jonswap" ? wave_spectrum::jonswap : wave_spectrum::pierson_moskowitz,
                                    it.need("Hs"), Tp, it.get("gamma", type=="jonswap" ? 3.3 : 1.0), dir,
                                    it.get("wmin", 0.5*wp), it.get("wmax", 3.0*wp), (int)it.get("n", 200),
                                    (unsigned)it.get("seed", 1));
            }
            else
                it.fail("unknown wave type '" + type + "'");
        }
    }
    // wind inflow: steady or uniform in time (wind steady / wind file) with IEC events and turbulence, or a
    // TurbSim full field (wind turbsim)
    if(!wind_items.empty())
    {
        const item *base = nullptr, *ts = nullptr, *iec = nullptr, *ev = nullptr, *turb = nullptr;
        for(const item &it : wind_items)
        {
            const std::string &type = it.t[1];
            const item **slot = type=="steady" || type=="file" ? &base : type=="turbsim" ? &ts : type=="iec" ? &iec
                                : type=="event" ? &ev : type=="turbulence" ? &turb : nullptr;
            if(!slot) it.fail("unknown wind type '" + type + "'");
            if(*slot) it.fail("wind " + type + " given twice (or wind steady and wind file)");
            *slot = &it;
        }
        iec_turbine tc;
        if(iec)
        {
            if(iec->t.size() < 4 || !tc.set_class(iec->t[2], iec->t[3]))
                iec->fail("wind iec <I|II|III> <A+|A|B|C> z_hub <z> D <D>");
            tc.z_hub = iec->need("z_hub");
            tc.D = iec->need("D");
        }
        if(ts)
        {
            if(base || ev || turb) ts->fail("wind turbsim is a full field: no wind steady, file, event or turbulence with it");
            auto g = std::make_shared<wind_grid>();
            std::string msg;
            if(!g->read_bts(path(ts->word(2)), msg)) ts->fail(msg);
            g->direction = ts->get("direction", 0.0)*deg;
            const int io = ts->find("origin");
            if(io >= 0) g->origin = Vec3(ts->num(io+1), ts->num(io+2), 0.0);
            g->xi0 = ts->get("x0", 0.0);
            wind = g;
        }
        else
        {
            if(!base) wind_items[0].fail("wind steady or wind file is needed (or wind turbsim)");
            const item &b = *base;
            const double zr = b.get("z_ref", iec ? tc.z_hub : -1.0);
            if(zr <= 0.0) b.fail("'z_ref <z>' is required (or wind iec ... z_hub)");
            auto w = std::make_shared<wind_field>(b.t[1]=="steady" ? b.need("U") : 0.0, zr, b.get("direction", 0.0)*deg);
            if(b.t[1]=="file")
            {
                std::vector<wind_field::series_row> rows;
                std::string msg;
                if(!wind_field::read_uniform(path(b.word(2)), rows, msg)) b.fail(msg);
                w->set_series(rows, b.get("L_ref", iec ? tc.D : 1.0));
            }
            if(b.has("shear") && b.has("log")) b.fail("shear (power law) or log, not both");
            if(b.has("shear")) w->set_power_law(b.need("shear"));
            if(b.has("log")) w->set_log_law(b.need("log"));
            w->set_veer(b.get("veer", 0.0)*deg);
            w->set_upflow(b.get("upflow", 0.0)*deg);
            w->set_ramp(b.get("ramp", 0.0));
            const int io = b.find("origin");
            if(io >= 0) w->set_origin(b.num(io+1), b.num(io+2));
            const double Vh = w->speed(0.0, iec ? tc.z_hub : zr, 0.0);
            if(ev)
            {
                if(!iec) ev->fail("wind event needs wind iec");
                const iec_event e = iec_event_from(ev->word(2));
                if(e==iec_event::none) ev->fail("unknown event '" + ev->word(2) + "' (eog, edc+, edc-, ecd+, ecd-, ews_v+, ews_v-, ews_h+, ews_h-)");
                w->set_event(e, tc, ev->need("t0"));
            }
            if(turb)
            {
                const item &u = *turb;
                if(!iec) u.fail("wind turbulence needs wind iec");
                if(u.word(2) != "kaimal") u.fail("turbulence model must be kaimal");
                const std::string tm = u.t.size() > 3 && (u.t[3]=="ntm" || u.t[3]=="etm" || u.t[3]=="ewm") ? u.t[3] : "ntm";
                const double s1 = tm=="ntm" ? tc.sigma_ntm(Vh) : tm=="etm" ? tc.sigma_etm(Vh) : tc.sigma_ewm(Vh);
                wind_grid::kaimal_settings ks;
                ks.ny = (int)u.get("ny", ks.ny); ks.nz = (int)u.get("nz", ks.nz);
                ks.width = u.get("width", 1.1*tc.D); ks.height = u.get("height", 1.1*tc.D);
                ks.T = u.get("T", run.t_end + 10.0); ks.dt = u.get("dt", ks.dt);
                ks.seed = (unsigned)u.get("seed", 1.0);
                ks.coherence_vw = u.has("coherence_vw");
                ks.scale_hub = u.has("scale_hub");
                if(ks.ny < 2 || ks.nz < 2 || ks.dt <= 0.0 || ks.T <= 0.0) u.fail("ny, nz >= 2, dt > 0 and T > 0 needed");
                if(tc.z_hub - 0.5*ks.height <= 0.0) u.fail("the grid reaches below z = 0 (height)");
                if(Vh <= 0.0) u.fail("the hub speed must be positive");
                auto g = std::make_shared<wind_grid>(wind_grid::kaimal(tc, Vh, s1, ks));
                if(u.has("write"))
                {
                    std::string msg;
                    if(!g->write_bts(path(u.word("write")), "CORAL IEC Kaimal " + tm, msg, w.get())) u.fail(msg);
                }
                w->set_turbulence(g);
            }
            wind = w;
        }
    }
    // the kinematics of Morison members and mooring lines: CORAL's waves, or the resolved fluid of REEF3D
    // sampled at their points (with CORAL's waves outside the fluid domain and before the first sample)
    std::shared_ptr<const wave_kinematics> kin = waves;
    if(fluid_kin_order >= 0)
    {
        fluid_kin = std::make_shared<sampled_kinematics>(waves);
        fluid_kin->order = fluid_kin_order;
        kin = fluid_kin;
    }
    // members (nodes are bodies: before the joints and loads), superelements, point masses, soil
    std::map<std::string, std::vector<int>> stations;     // frame member -> node bodies
    std::map<std::string, member_spec> frame_specs;
    std::map<std::string, superelement*> super_of;        // superelement member -> its superelement
    {
        auto parse_member = [&](const item &it)
        {
            member_spec ms;
            ms.name = it.word(1);
            ms.p1 = it.need_vec("p1");
            ms.p2 = it.need_vec("p2");
            ms.n = (int)it.get("n", 10);
            const std::string s1 = it.word("section");
            if(!sections.count(s1)) it.fail("unknown section '" + s1 + "'");
            ms.s1 = sections[s1];
            if(it.has("section2"))
            {
                const std::string s2 = it.word("section2");
                if(!sections.count(s2)) it.fail("unknown section '" + s2 + "'");
                ms.s2 = sections[s2];
                ms.tapered = true;
            }
            ms.orient = it.vec("orient", Vec3::UnitZ());
            ms.beta = it.get("beta", 0.0);
            if(it.has("end1")) ms.body1 = body_of(it, it.word("end1"));
            if(it.has("end2")) ms.body2 = body_of(it, it.word("end2"));
            return ms;
        };
        // members that belong to a superelement
        std::map<std::string, std::string> in_super;
        for(const item &it : items)
        {
            if(it.t[0] != "superelement") continue;
            const int i0 = it.find("members");
            if(i0 < 0) it.fail("'members' is required");
            for(int i=i0+1; i<(int)it.t.size() && it.t[i]!="modes" && it.t[i]!="zeta" && it.t[i]!="beta"; ++i)
            {
                if(in_super.count(it.t[i])) it.fail("member '" + it.t[i] + "' in two superelements");
                in_super[it.t[i]] = it.word(1);
            }
        }
        frame_builder fb(m);
        bool any = false;
        std::map<std::string, member_spec> super_members;
        for(const item &it : items)
        {
            if(it.t[0] != "member") continue;
            const member_spec ms = parse_member(it);
            if(stations.count(ms.name) || super_members.count(ms.name)) it.fail("member '" + ms.name + "' defined twice");
            if(in_super.count(ms.name))
            {
                super_members[ms.name] = ms;
                continue;
            }
            stations[ms.name] = fb.add_member(ms);
            frame_specs[ms.name] = ms;
            any = true;
        }
        if(any) fb.finish();
        for(const item &it : items)
        {
            if(it.t[0] != "superelement") continue;
            std::vector<member_spec> list;
            const int i0 = it.find("members");
            for(int i=i0+1; i<(int)it.t.size() && it.t[i]!="modes" && it.t[i]!="zeta" && it.t[i]!="beta"; ++i)
            {
                if(!super_members.count(it.t[i])) it.fail("unknown member '" + it.t[i] + "'");
                list.push_back(super_members[it.t[i]]);
            }
            try
            {
                auto se = std::make_unique<superelement>(m, it.word(1), list, (int)it.get("modes", 10), it.get("zeta", 0.0),
                                                         it.get("beta", 0.0), gravity_on ? m.gravity : Vec3::Zero());
                for(const auto &ms : list) super_of[ms.name] = se.get();
                m.add_load(std::move(se));
            }
            catch(const std::exception &e) {it.fail(e.what());}
        }
        for(const item &it : items)
        {
            if(it.t[0] != "point_mass") continue;
            const int b = body_of(it, it.word(1));
            if(b < 0) it.fail("point mass on the ground");
            Mat3 J = Mat3::Zero();
            const int ii = it.find("inertia");
            if(ii >= 0) J = it.vec(ii+1).asDiagonal();
            if(it.has("at"))
                m.add_point_mass(b, it.need_vec("at"), it.need("mass"), J);
            else
            {
                m.bodies[b].mass += it.need("mass");
                m.bodies[b].inertia += J;
            }
        }
        // soil profiles and pile springs
        std::map<std::string, soil_profile> soils;
        for(const item &it : items)
        {
            if(it.t[0] != "soil") continue;
            soil_profile sp;
            sp.name = it.word(1);
            sp.mudline = it.need("mudline");
            soils[sp.name] = sp;
        }
        for(const item &it : items)
        {
            if(it.t[0] != "soil_layer") continue;
            const std::string sn = it.word(1);
            if(!soils.count(sn)) it.fail("unknown soil '" + sn + "'");
            soil_layer l;
            const std::string type = it.word(2);
            l.top = it.need("top");
            l.bottom = it.need("bottom");
            l.gamma = it.get("gamma", l.gamma);
            l.cyclic = it.has("cyclic");
            if(type=="sand")
            {
                l.type = soil_layer::kind::sand;
                l.phi = it.need("phi");
                l.k = it.get("k", 0.0);
                l.delta = it.get("delta", -1.0);
                l.K = it.get("K", l.K);
                l.Nq = it.get("Nq", l.Nq);
            }
            else if(type=="clay")
            {
                l.type = soil_layer::kind::clay;
                l.su_top = it.need("su");
                l.su_bottom = it.get("su_bottom", -1.0);
                l.eps50 = it.get("eps50", l.eps50);
                l.J = it.get("J", l.J);
                l.residual = it.get("residual", l.residual);
            }
            else if(type=="linear")
            {
                l.type = soil_layer::kind::linear;
                l.kh = it.need("kh");
                l.kv = it.get("kv", 0.0);
                l.kq = it.get("kq", 0.0);
            }
            else
                it.fail("unknown soil layer type '" + type + "'");
            if(l.bottom <= l.top) it.fail("bottom must lie below top (depths below the mudline)");
            soils[sn].layers.push_back(l);
        }
        for(const item &it : items)
        {
            if(it.t[0] != "pile") continue;
            const std::string mem = it.word(1), sn = it.word("soil");
            if(!stations.count(mem)) it.fail("unknown frame member '" + mem + "'");
            if(!soils.count(sn)) it.fail("unknown soil '" + sn + "'");
            const bool all = !it.has("py") && !it.has("tz") && !it.has("qz");
            std::vector<int> nodes;
            for(int b : stations[mem]) if(b >= 0) nodes.push_back(b);
            try
            {
                add_pile_springs(m, soils[sn], nodes, it.need("D"), all || it.has("py"), all || it.has("tz"),
                                 all || it.has("qz"), it.has("linear"), "pile:" + mem);
            }
            catch(const std::exception &e) {it.fail(e.what());}
        }
    }

    if(gravity_on && m.nbody() > 0)
        m.add_load(std::make_unique<gravity_load>());

    // pass 2: joints
    for(const item &it : items)
    {
        if(it.t[0] != "joint") continue;
        const std::string name = it.word(1), type = it.word(2);
        joint_type jt;
        if(type=="fixed") jt = joint_type::fixed;
        else if(type=="spherical") jt = joint_type::spherical;
        else if(type=="revolute") jt = joint_type::revolute;
        else if(type=="prismatic") jt = joint_type::prismatic;
        else if(type=="cylindrical") jt = joint_type::cylindrical;
        else if(type=="planar") jt = joint_type::planar;
        else it.fail("unknown joint type '" + type + "'");
        const int a = body_of(it, it.word(3)), b = body_of(it, it.word(4));
        if(a==b) it.fail("joint between a body and itself");
        if(joints.count(name)) it.fail("joint '" + name + "' defined twice");
        joint &j = static_cast<joint&>(m.add_constraint(std::make_unique<joint>(m, jt, a, b, it.need_vec("point"), it.vec("axis", Vec3::UnitZ()), name)));
        joints[name] = &j;
    }

    // bathymetry of the dynamic mooring lines: plane, grid file, or the resolved fluid's bed (with
    // fluid_kinematics)
    std::shared_ptr<const seabed_surface> bathy;
    for(const item &it : items)
        if(it.t[0]=="bathymetry")
        {
            const std::string type = it.word(1);
            if(type=="plane")
                bathy = std::make_shared<seabed_plane>(it.num(2), it.t.size() > 3 ? it.num(3) : 0.0, it.t.size() > 4 ? it.num(4) : 0.0);
            else if(type=="grid")
                bathy = std::make_shared<seabed_grid>(seabed_grid::read(path(it.word(2))));
            else if(type=="fluid")
            {
                if(!fluid_kin) it.fail("bathymetry fluid needs fluid_kinematics");
                fluid_bed = true;
            }
            else
                it.fail("bathymetry plane <z0> [<sx> <sy>] | grid <file> | fluid");
        }
    // line types: properties of mooring lines for 'type' and 'section'
    std::map<std::string, mooring_props> line_types;
    auto read_props = [&](const item &it, mooring_props p, bool need)
    {
        if(need || it.has("d")) p.d = need ? it.need("d") : it.get("d", p.d);
        if(need || it.has("mass")) p.mass = need ? it.need("mass") : it.get("mass", p.mass);
        if(need || it.has("EA")) p.EA = need ? it.need("EA") : it.get("EA", p.EA);
        p.BA = it.get("BA", p.BA); p.CB = it.get("CB", p.CB);
        p.Cdn = it.get("Cdn", p.Cdn); p.Cdt = it.get("Cdt", p.Cdt); p.Can = it.get("Can", p.Can); p.Cat = it.get("Cat", p.Cat);
        p.kb = it.get("kb", p.kb); p.cb = it.get("cb", p.cb); p.mu = it.get("mu", p.mu);
        p.mu_axial = it.get("mu_axial", p.mu_axial); p.v_eps = it.get("v_eps", p.v_eps);
        p.mu_static = it.get("mu_static", p.mu_static); p.ks = it.get("ks", p.ks); p.cs = it.get("cs", p.cs);
        if(!(p.mu_static > 0.0)) it.fail("mu_static (static / kinetic friction) must be > 0");
        p.EI = it.get("EI", p.EI);
        if(it.has("BA_slack")) p.BA_slack = true;
        if(p.v_eps < 0.0) it.fail("v_eps must be >= 0 (0: stick-slip friction)");
        if(p.d <= 0.0 || p.EA <= 0.0 || p.mass < 0.0 || p.EI < 0.0) it.fail("d and EA must be > 0, mass and EI >= 0");
        p.rho = rho;
        return p;
    };
    for(const item &it : items)
        if(it.t[0]=="line_type")
        {
            const std::string n = it.word(1);
            if(line_types.count(n)) it.fail("line type '" + n + "' is defined twice");
            line_types[n] = read_props(it, mooring_props(), true);
        }
    auto line_type = [&](const item &it, const std::string &n) -> const mooring_props&
    {
        const auto q = line_types.find(n);
        if(q == line_types.end()) it.fail("unknown line type '" + n + "'");
        return q->second;
    };
    // the sections of a line: 'section <type> <L> <nseg>' (in order from a to b)
    struct line_section {mooring_props p; int nseg;};
    auto sections_of = [&](const item &it)
    {
        std::vector<line_section> v;
        for(int i=it.find("section"); i >= 0; i=it.find("section", i+1))
        {
            line_section q{line_type(it, it.word(i+1)), (int)it.num(i+3)};
            q.p.L = it.num(i+2);
            if(q.p.L <= 0.0 || q.nseg < 1) it.fail("section: length > 0 and at least one segment");
            v.push_back(q);
        }
        return v;
    };

    // mooring systems: the free points (point) and the dynamic lines that end on them, grouped by their
    // connections (a line of several sections is a group with its junctions as points); a group is built
    // where its first line stands in the file
    struct moor_group {std::vector<const item*> lines; std::vector<const item*> points; bool built = false;};
    std::vector<moor_group> groups;
    std::map<const item*, int> group_of_line;
    std::map<std::string, const item*> point_item;
    {
        std::vector<const item*> pts;
        for(const item &it : items)
            if(it.t[0]=="point")
            {
                const std::string n = it.word(1);
                if(point_item.count(n)) it.fail("point '" + n + "' is defined twice");
                if(n=="ground" || m.find_body(n) >= 0) it.fail("point '" + n + "': a body has that name");
                point_item[n] = &it;
                pts.push_back(&it);
            }
        // union-find over the points; the lines join the points they end on
        std::map<std::string, std::string> parent;
        std::function<std::string(const std::string&)> root = [&](const std::string &x) -> std::string
        {
            return parent[x]==x ? x : (parent[x] = root(parent[x]));
        };
        for(const auto &q : point_item) parent[q.first] = q.first;
        std::vector<std::pair<const item*, std::vector<std::string>>> dl;
        for(const item &it : items)
            if(it.t[0]=="mooring" && it.t.size() > 2 && it.t[2]=="dynamic")
            {
                std::vector<std::string> pe;
                for(const char *key : {"a", "b"})
                {
                    const int i = it.find(key);
                    if(i >= 0 && i+1 < (int)it.t.size() && point_item.count(it.t[i+1]))
                        pe.push_back(it.t[i+1]);
                }
                if(sections_of(it).size() > 1)
                {
                    // the junctions of the sections: a group of its own (joined with the points at the ends)
                    const std::string j = "\x01" + it.word(1);
                    parent[j] = j;
                    pe.push_back(j);
                }
                for(size_t q=1; q<pe.size(); ++q) parent[root(pe[0])] = root(pe[q]);
                if(!pe.empty()) dl.push_back({&it, pe});
            }
        std::map<std::string, int> gid;
        for(const auto &x : dl)
        {
            const std::string r = root(x.second[0]);
            if(!gid.count(r)) {gid[r] = (int)groups.size(); groups.emplace_back();}
            groups[gid[r]].lines.push_back(x.first);
            group_of_line[x.first] = gid[r];
        }
        for(const item *it : pts)
        {
            const std::string r = root(it->word(1));
            if(!gid.count(r)) it->fail("point '" + it->word(1) + "' is not an end of a dynamic mooring line");
            groups[gid[r]].points.push_back(it);
        }
    }
    // the properties of a line: its own keys, over those of 'type <line type>' if given
    auto line_props = [&](const item &it)
    {
        const bool typed = it.has("type");
        mooring_props p = read_props(it, typed ? line_type(it, it.word("type")) : mooring_props(), !typed);
        p.L = it.need("L");
        return p;
    };
    // an end 'a <body> <x y z>' or 'a <point>' (point: index in the group's points, else -1)
    auto line_end = [&](const item &it, const char *key, const std::map<std::string, int> &pidx)
    {
        const int i = it.find(key);
        if(i < 0) it.fail("'a <body> <x y z>' (or 'a <point>') and 'b <body> <x y z>' are required");
        mooring_end e;
        const std::string n = it.word(i+1);
        const auto q = pidx.find(n);
        if(q != pidx.end())
            e.point = q->second;
        else
        {
            e.body = body_of(it, n);
            e.p = it.vec(i+2);
        }
        return e;
    };
    // bending: clamp_a [<dx dy dz>], clamp_b [...]: the end clamped to its body along the direction (default: the
    // start tangent of the line there)
    auto clamps = [&](const item &it, mooring_system &sys, int ka, int kb)
    {
        for(int e=0; e<2; ++e)
        {
            const int i = it.find(e==0 ? "clamp_a" : "clamp_b");
            if(i < 0)
                continue;
            Vec3 d = Vec3::Zero();
            auto isnum = [&](int j) {if(j >= (int)it.t.size()) return false; try {size_t n; std::stod(it.t[j], &n); return n==it.t[j].size();} catch(...) {return false;}};
            if(isnum(i+1) && isnum(i+2) && isnum(i+3))
                d = it.vec(i+1);
            try {sys.set_clamp(e==0 ? ka : kb, e, m, d);}
            catch(const std::exception &x) {it.fail(x.what());}
        }
    };
    // failure of a dynamic line: break_time <t>, break_tension <T>, break_end <a|b>
    auto line_break = [&](const item &it, bool &any)
    {
        mooring_break b;
        b.t = it.get("break_time", b.t);
        b.T = it.get("break_tension", b.T);
        if(it.has("break_end"))
        {
            const std::string e = it.word("break_end");
            if(e=="a") b.end = 0;
            else if(e=="b") b.end = 1;
            else it.fail("break_end must be a or b");
        }
        any = it.has("break_time") || it.has("break_tension");
        if(any && !(b.T > 0.0)) it.fail("break_tension must be > 0");
        return b;
    };
    auto build_group = [&](moor_group &g)
    {
        const item &first = *g.lines[0];
        auto sys = std::make_unique<mooring_system>(m, kin, first.vec("current", Vec3::Zero()), first.word(1));
        if(bathy || fluid_bed)
            sys->set_seabed(bathy, fluid_bed);
        std::map<std::string, int> pidx;
        // the seabed of the group: the lowest ground end of its lines
        double zsys = -1.0e30;
        bool any = false;
        for(const item *li : g.lines)
            for(const char *key : {"a", "b"})
            {
                const int i = li->find(key);
                const std::string n = li->word(i+1);
                if(!point_item.count(n) && body_of(*li, n) < 0)
                {
                    const double z = li->num(i+4);
                    zsys = any ? std::min(zsys, z) : z;
                    any = true;
                }
            }
        for(const item *pi : g.points)
        {
            mooring_point_props pp;
            pp.mass = pi->get("mass", 0.0);
            pp.volume = pi->get("volume", 0.0);
            pp.CdA = pi->get("CdA", 0.0);
            pp.Ca = pi->get("Ca", 0.0);
            pp.radius = pi->get("radius", -1.0);
            pp.kb = pi->get("kb", pp.kb);
            pp.cb = pi->get("cb", pp.cb);
            pp.mu = pi->get("mu", 0.0);
            pp.v_eps = pi->get("v_eps", pp.v_eps);
            pp.mu_static = pi->get("mu_static", pp.mu_static);
            pp.ks = pi->get("ks", pp.ks);
            pp.cs = pi->get("cs", pp.cs);
            if(pp.v_eps < 0.0 || !(pp.mu_static > 0.0)) pi->fail("v_eps must be >= 0 (0: stick-slip friction), mu_static > 0");
            if(pp.mass < 0.0 || pp.volume < 0.0 || pp.CdA < 0.0) pi->fail("mass, volume and CdA must not be negative");
            const double zb = pi->has("no_seabed") ? -1.0e30 : pi->get("seabed", zsys);
            pidx[pi->word(1)] = sys->add_point(pi->word(1), pi->vec(2), pp, zb);
            if(pi->has("seabed")) sys->set_point_flat(pidx[pi->word(1)]);
        }
        double dt = 0.0, settle = 0.0;
        bool exact = false, nostatic = false;
        for(const item *li : g.lines)
        {
            const item &it = *li;
            const mooring_end ea = line_end(it, "a", pidx), eb = line_end(it, "b", pidx);
            const bool ground = (ea.body < 0 && ea.point < 0) || (eb.body < 0 && eb.point < 0);
            const double zdef = ground ? default_seabed(ea, eb) : zsys;
            const double zbed = it.has("no_seabed") ? -1.0e30 : it.get("seabed", zdef);
            if((it.vec("current", Vec3::Zero()) - first.vec("current", Vec3::Zero())).norm() > 0.0)
                it.fail("the lines of a mooring system take one current (that of " + first.word(1) + ")");
            const std::vector<line_section> sec = sections_of(it);
            if(sec.size() > 1)
            {
                // junction points between the sections (no mass, volume or drag of their own), the sections
                // as lines '<name>.1' .. '<name>.n', the junctions '<name>.p1' ..
                auto pos = [&](const mooring_end &e) {return e.point >= 0 ? sys->point_position(e.point) : e.p;};
                std::vector<mooring_props> sp;
                for(const auto &q : sec) sp.push_back(q.p);
                const std::vector<Vec3> J = mooring_junctions(pos(ea), pos(eb), sp, zbed, m.gravity.norm());
                std::vector<mooring_end> ends(sec.size() + 1);
                ends[0] = ea;
                ends[sec.size()] = eb;
                for(size_t q=0; q<J.size(); ++q)
                {
                    ends[q+1].point = sys->add_point(it.word(1) + ".p" + std::to_string(q+1), J[q], mooring_point_props(), zbed);
                    if(it.has("seabed")) sys->set_point_flat(ends[q+1].point);
                }
                int first = -1, last = -1;
                for(size_t q=0; q<sec.size(); ++q)
                {
                    last = sys->add_line(m, it.word(1) + "." + std::to_string(q+1), ends[q], ends[q+1], sec[q].p, sec[q].nseg, zbed);
                    if(it.has("seabed")) sys->set_line_flat(last);
                    if(q==0) first = last;
                }
                clamps(it, *sys, first, last);
                // a failure of the line: of the section at that end (the tension: at either end)
                bool any = false;
                const mooring_break b = line_break(it, any);
                if(any)
                {
                    mooring_break ba = b, bb = b;
                    ba.end = 0;
                    bb.end = 1;
                    if(b.end != 1) {if(b.end < 0) ba.t = std::numeric_limits<double>::infinity(); sys->set_break(first, ba);}
                    if(b.end != 0) sys->set_break(last, bb);
                }
            }
            else
            {
                const int q = sys->add_line(m, it.word(1), ea, eb, sec.size()==1 ? sec[0].p : line_props(it),
                                            sec.size()==1 ? sec[0].nseg : (int)it.get("nseg", 20), zbed);
                if(it.has("seabed")) sys->set_line_flat(q);
                clamps(it, *sys, q, q);
                bool any = false;
                const mooring_break b = line_break(it, any);
                if(any) sys->set_break(q, b);
            }
            const double d = it.get("dt_line", 0.0);
            if(d > 0.0) dt = dt > 0.0 ? std::min(dt, d) : d;
            exact = exact || it.has("exact_kinematics");
            nostatic = nostatic || it.has("no_static");
            if(it.has("settle")) settle = std::max(settle, it.need("settle"));
        }
        sys->set_dt_line(dt);
        sys->set_exact_kinematics(exact);
        sys->set_static_start(!nostatic);
        if(settle > 0.0)
            settles.push_back({sys.get(), settle});
        m.add_load(std::move(sys));
        g.built = true;
    };

    // pass 3: loads in file order
    for(const item &it : items)
    {
        const std::string &k = it.t[0];
        if(k=="hydro")
        {
            const std::string name = it.word(1);
            auto db = std::make_shared<hydro_db>();
            const int ib = it.find("bodies");
            if(ib < 0) it.fail("'bodies' is required");
            std::vector<int> map;
            for(int i=ib+1; i<(int)it.t.size(); ++i)
            {
                const int b = m.find_body(it.t[i]);
                if(b < 0) break;
                map.push_back(b);
            }
            if(map.empty()) it.fail("no known body after 'bodies'");
            if(it.has("file"))
                db->read(path(it.word("file")));
            else if(it.has("wamit"))
                db->read_wamit(path(it.word("wamit")), (int)map.size(), it.need("ulen"), rho, std::fabs(m.gravity(2)), waves ? waves->depth() : -1.0);
            else
                it.fail("'file' or 'wamit' is required");
            if(db->nbody() != (int)map.size())
                it.fail("the database has " + std::to_string(db->nbody()) + " bodies, 'bodies' lists " + std::to_string(map.size()));
            for(int i=it.find("volume"); i>=0; i=it.find("volume", i+1))
                db->body.at((int)it.num(i+1) - 1).volume = it.num(i+2);
            for(int i=it.find("ref"); i>=0; i=it.find("ref", i+1))
                db->body.at((int)it.num(i+1) - 1).ref = it.vec(i+2);
            hydro_settings s;
            s.state_space = it.has("state_space");
            s.radiation = !it.has("no_radiation");
            s.excitation = !it.has("no_excitation");
            s.hydrostatics = !it.has("no_hydrostatics");
            s.buoyancy = !it.has("no_buoyancy");
            const int im = it.find("memory");
            if(im >= 0) {s.t_memory = it.num(im+1); s.dt_memory = it.num(im+2);}
            const int nd = db->ndof();
            for(int i=it.find("damping"); i>=0; i=it.find("damping", i+1))
            {
                if(s.linear_damping.size()==0) s.linear_damping = MatX::Zero(nd, nd);
                const int d = (int)it.num(i+1) - 1;
                if(d < 0 || d >= nd) it.fail("damping mode out of range");
                s.linear_damping(d,d) = it.num(i+2);
            }
            for(int i=it.find("drag"); i>=0; i=it.find("drag", i+1))
            {
                if(s.quadratic_drag.size()==0) s.quadratic_drag = MatX::Zero(nd, nd);
                const int d = (int)it.num(i+1) - 1;
                if(d < 0 || d >= nd) it.fail("drag mode out of range");
                s.quadratic_drag(d,d) = it.num(i+2);
            }
            m.add_load(std::make_unique<hydro_load>(m, db, map, waves, s, name));
        }
        else if(k=="pto")
        {
            const std::string name = it.word(1), jn = it.word("joint");
            if(!joints.count(jn)) it.fail("unknown joint '" + jn + "'");
            pto_settings s;
            s.k = it.get("k", 0.0);
            s.c = it.get("c", 0.0);
            s.q0 = it.get("q0", 0.0);
            s.F0 = it.get("F0", 0.0);
            s.coulomb = it.get("coulomb", 0.0);
            s.v_eps = it.get("v_eps", s.v_eps);
            s.c_loss = it.get("c_loss", 0.0);
            s.coulomb_loss = it.get("coulomb_loss", 0.0);
            if(s.c_loss < 0.0 || s.coulomb_loss < 0.0) it.fail("c_loss and coulomb_loss must be >= 0");
            s.F_max = it.get("F_max", 0.0);
            s.P_max = it.get("P_max", 0.0);
            const int is = it.find("stops");
            if(is >= 0)
            {
                s.q_min = it.num(is+1);
                s.q_max = it.num(is+2);
                s.k_stop = it.num(is+3);
                if(is+4 < (int)it.t.size()) {try {s.c_stop = std::stod(it.t[is+4]);} catch(...) {}}
            }
            const int il = it.find("latching");
            if(il >= 0)
            {
                s.control = pto_control::latching;
                s.t_hold = it.num(il+1);
                s.c_latch = it.num(il+2);
                if(il+3 < (int)it.t.size()) {try {s.k_latch = std::stod(it.t[il+3]);} catch(...) {}}
            }
            const int idc = it.find("declutching");
            if(idc >= 0) {s.control = pto_control::declutching; s.t_hold = it.num(idc+1);}
            m.add_load(std::make_unique<pto_load>(*joints[jn], it.has("rotational"), s, name));
        }
        else if(k=="drive")
        {
            // prescribed motion of a joint coordinate (forced oscillation)
            const std::string name = it.word(1), jn = it.word("joint");
            if(!joints.count(jn)) it.fail("unknown joint '" + jn + "'");
            drive_settings s;
            const int is = it.find("sine");
            if(is < 0) it.fail("'sine <amplitude> <period>' is required");
            s.amp = it.num(is+1);
            s.period = it.num(is+2);
            if(s.period <= 0.0) it.fail("the period must be positive");
            s.q0 = it.get("q0", 0.0);
            s.phase = it.get("phase", 0.0);
            s.ramp = it.get("ramp", 0.0);
            s.k = it.get("k", s.k);
            s.c = it.get("c", s.c);
            m.add_load(std::make_unique<joint_drive>(*joints[jn], it.has("rotational"), s, name));
        }
        else if(k=="jforce")
        {
            // prescribed force along a joint axis from a table t F (system identification)
            const std::string name = it.word(1), jn = it.word("joint");
            if(!joints.count(jn)) it.fail("unknown joint '" + jn + "'");
            jforce_settings s;
            std::string msg;
            if(!joint_force::read_table(path(it.word("file")), s.t, s.F, msg)) it.fail(msg);
            s.scale = it.get("scale", 1.0);
            s.ramp = it.get("ramp", 0.0);
            s.periodic = it.has("periodic");
            m.add_load(std::make_unique<joint_force>(*joints[jn], it.has("rotational"), s, name));
        }
        else if(k=="spring")
        {
            const int ia = it.find("a"), ib = it.find("b");
            if(ia < 0 || ib < 0) it.fail("'a <body> <x y z>' and 'b <body> <x y z>' are required");
            const int a = body_of(it, it.word(ia+1)), b = body_of(it, it.word(ib+1));
            m.add_load(std::make_unique<spring_damper>(m, a, it.vec(ia+2), b, it.vec(ib+2), it.need("k"), it.get("c", 0.0),
                                                       it.get("L0", -1.0), it.get("T0", 0.0), it.has("tension_only"), it.word(1)));
        }
        else if(k=="linear6")
        {
            const int b = body_of(it, it.word("body"));
            Mat6 K = Mat6::Zero(), C = Mat6::Zero();
            const int ik = it.find("k"), ic = it.find("c"), iK = it.find("K"), iC = it.find("C");
            for(int d=0; d<6; ++d)
            {
                if(ik >= 0) K(d,d) = it.num(ik+1+d);
                if(ic >= 0) C(d,d) = it.num(ic+1+d);
            }
            // full matrices, row by row
            for(int d=0; d<36; ++d)
            {
                if(iK >= 0) K(d/6, d%6) = it.num(iK+1+d);
                if(iC >= 0) C(d/6, d%6) = it.num(iC+1+d);
            }
            m.add_load(std::make_unique<linear6>(m, b, it.need_vec("point"), K, C, Vec6::Zero(), it.word(1)));
        }
        else if(k=="drag")
        {
            // quadratic damping of a body (coral_body_drag.h)
            const int b = body_of(it, it.word("body"));
            if(b < 0) it.fail("the ground has no drag");
            const Vec3 P = it.vec("point", m.body(b).x);
            const double r = it.get("rho", rho);
            const int iq = it.find("Q");
            if(iq >= 0)
            {
                Vec6 Q;
                for(int d=0; d<6; ++d)
                    Q(d) = it.num(iq+1+d);
                m.add_load(std::make_unique<body_drag>(m, b, P, Q, it.word(1)));
            }
            else
            {
                const Vec3 e = it.vec("axis", Vec3(0.0, 0.0, 1.0));
                if(e.norm() <= 0.0) it.fail("axis must not be zero");
                m.add_load(std::make_unique<body_drag>(m, b, P, e, 0.5*r*it.need("Cd")*it.need("area"), it.word(1)));
            }
        }
        else if(k=="mooring")
        {
            const std::string name = it.word(1), type = it.word(2);
            if(group_of_line.count(&it))
            {
                moor_group &g = groups[group_of_line[&it]];
                if(!g.built)
                    build_group(g);
                continue;
            }
            const int ia = it.find("a"), ib = it.find("b");
            if(ia < 0 || ib < 0) it.fail("'a <body> <x y z>' and 'b <body> <x y z>' are required");
            mooring_end ea, eb;
            ea.body = body_of(it, it.word(ia+1)); ea.p = it.vec(ia+2);
            eb.body = body_of(it, it.word(ib+1)); eb.p = it.vec(ib+2);
            const std::vector<line_section> sec = sections_of(it);
            if(sec.size() > 1) it.fail("a quasi-static line has one section (several: a dynamic line)");
            const mooring_props p = sec.size()==1 ? sec[0].p : line_props(it);
            // the seabed: at the ground end (the lower of two), none between two bodies (shared lines)
            const double zbed = it.has("no_seabed") ? -1.0e30 : it.get("seabed", default_seabed(ea, eb));
            if(type=="quasi_static")
                m.add_load(std::make_unique<mooring_qs>(m, ea, eb, p, name, zbed));
            else if(type=="dynamic")
            {
                auto l = std::make_unique<mooring_line>(m, ea, eb, p, sec.size()==1 ? sec[0].nseg : (int)it.get("nseg", 20), kin, it.vec("current", Vec3::Zero()),
                                                        zbed, it.get("dt_line", 0.0), name);
                l->set_exact_kinematics(it.has("exact_kinematics"));
                l->set_static_start(!it.has("no_static"));
                if((bathy || fluid_bed) && !it.has("seabed"))
                    l->set_seabed(bathy, fluid_bed);
                clamps(it, *l, 0, 0);
                bool any = false;
                const mooring_break b = line_break(it, any);
                if(any) l->set_break(0, b);
                if(it.has("settle"))
                    settles.push_back({l.get(), it.need("settle")});
                m.add_load(std::move(l));
            }
            else
                it.fail("mooring type must be quasi_static or dynamic");
        }
        else if(k=="force")
        {
            const int b = body_of(it, it.word("body"));
            auto pf = std::make_unique<point_force>(m, b, it.vec("point", m.body(b).x), it.vec("F", Vec3::Zero()),
                                                    it.vec("M", Vec3::Zero()), it.get("ramp", 0.0), it.word(1));
            pf->until = it.get("until", pf->until);
            m.add_load(std::move(pf));
        }
        else if(k=="morison" && it.has("member"))
        {
            // distributed Morison loads on a frame member or a superelement member
            const std::string mem = it.word("member");
            member_morison_coef c;
            c.Cd = it.get("Cd", 1.0);
            c.Ca = it.get("Ca", 1.0);
            c.Cm = it.get("Cm", 1.0 + c.Ca);
            c.n = (int)it.get("n", 2);
            c.buoyancy = !it.has("no_buoyancy");
            std::shared_ptr<const point_set> pts;
            std::vector<std::pair<int,int>> elems;
            const double g = std::fabs(m.gravity(2)) > 0.0 ? std::fabs(m.gravity(2)) : 9.81;
            if(stations.count(mem))
            {
                const member_spec &ms = frame_specs[mem];
                c.D = it.get("D", ms.s1.D);
                auto bp = std::make_shared<body_points>();
                const std::vector<int> &st = stations[mem];
                const int n = (int)st.size() - 1;
                auto on_body = [&](int i) {return (i==0 && ms.body1 > -2) || (i==n && ms.body2 > -2);};
                for(int i=0; i<=n; ++i)
                {
                    bp->add(m, st[i], ms.p1 + (double)i/n*(ms.p2 - ms.p1));
                    if(i > 0) elems.emplace_back(i-1, i);
                }
                pts = bp;
                // flooded member: the inner water below the still water level as node mass
                if(it.has("flooded"))
                {
                    const double Din = it.need("flooded"), Ai = 0.25*M_PI*Din*Din;
                    for(int i=0; i<n; ++i)
                    {
                        const Vec3 pa = ms.p1 + (double)i/n*(ms.p2 - ms.p1), pb = ms.p1 + (double)(i+1)/n*(ms.p2 - ms.p1);
                        double f = 0.0;
                        if(pa(2) <= 0.0 && pb(2) <= 0.0) f = 1.0;
                        else if(pa(2) < 0.0 || pb(2) < 0.0) f = -std::min(pa(2), pb(2))/std::fabs(pb(2) - pa(2));
                        const double mw = 0.5*rho*Ai*(pb - pa).norm()*f;
                        for(int e=0; e<2; ++e)
                        {
                            const int b = st[i+e];
                            if(b < 0 || mw <= 0.0) continue;
                            if(on_body(i+e)) m.add_point_mass(b, e ? pb : pa, mw);
                            else m.bodies[b].mass += mw;
                        }
                    }
                }
            }
            else if(super_of.count(mem))
            {
                superelement *se = super_of[mem];
                if(it.has("flooded")) it.fail("flooded superelement members: put the inner water into the section mass");
                c.D = it.get("D", se->member_diameter(mem));
                pts = std::shared_ptr<const point_set>(std::shared_ptr<const point_set>(), se);     // owned by the model
                const std::vector<int> nd = se->member_nodes(mem);
                for(size_t i=1; i<nd.size(); ++i) elems.emplace_back(nd[i-1], nd[i]);
            }
            else
                it.fail("unknown member '" + mem + "'");
            if(c.D <= 0.0) it.fail("'D' is required (the section has no diameter)");
            m.add_load(std::make_unique<member_morison>(pts, elems, c, kin, rho, g, it.vec("current", Vec3::Zero()), it.word(1)));
        }
        else if(k=="strip_member")
        {
            // strip theory with nonlinear hydrostatics for a rigid floater (coral_strip_hydro.h)
            strip_member sm;
            sm.p1 = it.need_vec("p1");
            sm.p2 = it.need_vec("p2");
            sm.D = it.need("D");
            sm.Cd = it.get("Cd", 1.0);
            sm.Ca = it.get("Ca", 1.0);
            sm.Cm = it.get("Cm", 1.0 + sm.Ca);
            sm.n = (int)it.get("n", 20);
            sm.end1 = it.has("end1");
            sm.end2 = it.has("end2");
            sm.Ca_end = it.get("Ca_end", 1.0);
            sm.Cd_end = it.get("Cd_end", 1.0);
            sm.buoyancy = !it.has("no_buoyancy");
            const int b = body_of(it, it.word("body"));
            if(b < 0) it.fail("strip_member needs a body");
            const double g = std::fabs(m.gravity(2)) > 0.0 ? std::fabs(m.gravity(2)) : 9.81;
            auto sh = std::make_unique<strip_hydro>(m, b, std::vector<strip_member>{sm}, kin, rho, g, it.vec("current", Vec3::Zero()), it.word(1));
            sh->exact_kinematics = it.has("exact_kinematics");
            m.add_load(std::move(sh));
        }
        else if(k=="morison")
        {
            morison_member mb;
            mb.p1 = it.need_vec("p1");
            mb.p2 = it.need_vec("p2");
            mb.D = it.need("D");
            mb.Cd = it.get("Cd", 1.0);
            mb.Ca = it.get("Ca", 1.0);
            mb.Cm = it.get("Cm", 1.0 + mb.Ca);
            mb.n = (int)it.get("n", 20);
            const int b = body_of(it, it.word("body"));
            m.add_load(std::make_unique<morison_load>(m, b, std::vector<morison_member>{mb}, kin, rho, it.vec("current", Vec3::Zero()), it.word(1)));
        }
        else if(k=="fluid_body")
        {
            // body driven by the resolved fluid of REEF3D (coral/reef3d): one fluid load for all of them
            fluid_body_spec fb;
            fb.body = body_of(it, it.word(1));
            if(fb.body < 0) it.fail("the ground cannot be a fluid body");
            fb.n = (int)it.get("n", (double)(fluid_bodies.size() + fluid_members.size()));
            fb.own_mass = it.has("own_mass");
            fb.one_way = it.has("one_way");
            fb.excitation = it.has("excitation");
            if(fb.excitation && !fb.one_way) it.fail("excitation needs one_way");
            const int ia = it.find("added_mass");
            if(ia >= 0)
            {
                if(fb.one_way) it.fail("added_mass is for two-way bodies, not one_way");
                for(int d=0; d<6; ++d)
                {
                    fb.added_mass(d) = it.num(ia+1+d);
                    if(fb.added_mass(d) < 0.0) it.fail("added_mass must not be negative");
                }
            }
            for(const auto &o : fluid_bodies)
            {
                if(o.body == fb.body) it.fail("body '" + it.word(1) + "' is a fluid body twice");
                if(o.n == fb.n) it.fail("REEF3D body " + std::to_string(fb.n) + " mapped twice");
            }
            for(const auto &o : fluid_members)
                if(o.n == fb.n && !fb.one_way) it.fail("REEF3D body " + std::to_string(fb.n) + " is a fluid_member: as fluid_body only one_way");
            fluid_bodies.push_back(fb);
        }
        else if(k=="fluid_member")
        {
            // distributed one-way loads of a REEF3D 6DOF body on a frame or superelement member
            const std::string mem = it.word("member");
            std::shared_ptr<const point_set> pts;
            std::vector<std::pair<int,int>> elems;
            if(stations.count(mem))
            {
                const member_spec &ms = frame_specs[mem];
                auto bp = std::make_shared<body_points>();
                const std::vector<int> &st = stations[mem];
                const int n = (int)st.size() - 1;
                for(int i=0; i<=n; ++i)
                {
                    bp->add(m, st[i], ms.p1 + (double)i/n*(ms.p2 - ms.p1));
                    if(i > 0) elems.emplace_back(i-1, i);
                }
                pts = bp;
            }
            else if(super_of.count(mem))
            {
                superelement *se = super_of[mem];
                pts = std::shared_ptr<const point_set>(std::shared_ptr<const point_set>(), se);
                const std::vector<int> nd = se->member_nodes(mem);
                for(size_t i=1; i<nd.size(); ++i) elems.emplace_back(nd[i-1], nd[i]);
            }
            else
                it.fail("unknown member '" + mem + "'");
            if(elems.empty()) it.fail("member '" + mem + "' has no elements");
            Vec3 pa, pb, v, a;
            pts->point_state(m, elems.front().first, pa, v, a);
            pts->point_state(m, elems.back().second, pb, v, a);
            const int ns = (int)it.get("strips", (double)(2*elems.size()));
            auto l = std::make_unique<member_fluid_load>(m, pts, elems, pa, pb, ns, it.word(1));
            l->order = (int)it.get("order", 2.0);
            l->excitation = it.has("excitation");
            fluid_member_spec fm;
            fm.load = l.get();
            fm.n = (int)it.get("n", (double)(fluid_bodies.size() + fluid_members.size()));
            for(const auto &o : fluid_bodies)
                if(o.n == fm.n && !o.one_way) it.fail("REEF3D body " + std::to_string(fm.n) + " is a two-way fluid_body");
            for(const auto &o : fluid_members)
                if(o.n == fm.n) it.fail("REEF3D body " + std::to_string(fm.n) + " mapped twice");
            fluid_members.push_back(fm);
            m.add_load(std::move(l));
        }
        else if(k=="rotor")
        {
            // wind turbine rotor with rigid blades, drivetrain and generator (coral/aero/coral_rotor.h)
            if(!wind) it.fail("rotor needs a wind inflow (wind ...)");
            rotor_spec sp;
            sp.name = it.word(1);
            sp.B = (int)it.get("blades", 3.0);
            sp.hub = it.need_vec("hub");
            sp.axis = it.need_vec("axis");
            if(sp.axis.norm() <= 0.0) it.fail("axis must not be zero");
            sp.R_hub = it.need("R_hub");
            std::string msg;
            if(!read_aerodyn_blade(path(it.word("blade")), sp.R_hub, sp.stations, msg)) it.fail(msg);
            const int ia = it.find("airfoils");
            if(ia < 0) it.fail("'airfoils <n> <file 1> .. <file n>' is required");
            const int na = (int)it.num(ia+1);
            for(int q=0; q<na; ++q)
            {
                airfoil af;
                if(!airfoil::read(path(it.word(ia+2+q)), af, msg)) it.fail(msg);
                sp.airfoils.push_back(af);
            }
            if(it.has("structure") && !blade_structure::read(path(it.word("structure")), sp.structure, msg)) it.fail(msg);
            sp.hub_mass = it.get("hub_mass", 0.0);
            const int ij = it.find("hub_inertia");
            if(ij >= 0)
            {
                sp.hub_J = it.num(ij+1);
                double v;
                if(ij+2 < (int)it.t.size()) {try {size_t n; v = std::stod(it.t[ij+2], &n); if(n==it.t[ij+2].size()) sp.hub_Jt = v;} catch(...) {}}
            }
            sp.precone = it.get("precone", 0.0)*deg;
            sp.azimuth0 = it.get("azimuth", 0.0)*deg;
            sp.ccw = it.has("ccw");
            if(it.has("nacelle")) sp.nacelle = body_of(it, it.word("nacelle"));
            sp.omega0 = it.get("rpm", 0.0)*M_PI/30.0;
            sp.pitch = it.get("pitch", 0.0)*deg;
            if(it.has("pitch_file"))
            {
                if(!joint_force::read_table(path(it.word("pitch_file")), sp.pitch_t, sp.pitch_v, msg)) it.fail(msg);
                for(double &x : sp.pitch_v) x *= deg;
            }
            sp.dt.N = it.get("gearbox", 1.0);
            sp.J_gen = it.get("gen_inertia", 0.0);
            const int is = it.find("shaft");
            if(is >= 0) {sp.dt.k = it.num(is+1); sp.dt.c = it.num(is+2);}
            sp.dt.eta = it.get("efficiency", 1.0);
            const int ig = it.find("generator");
            if(ig >= 0)
            {
                const std::string g = it.word(ig+1);
                if(g=="none") sp.dt.law = generator_law::none;
                else if(g=="constant") {sp.dt.law = generator_law::constant; sp.dt.T = it.num(ig+2);}
                else if(g=="slip") {sp.dt.law = generator_law::slip; sp.dt.slope = it.num(ig+2); sp.dt.omega_sync = it.num(ig+3)*M_PI/30.0;}
                else if(g=="k_omega2")
                {
                    sp.dt.law = generator_law::k_omega2;
                    sp.dt.K = it.num(ig+2);
                    double v;
                    if(ig+3 < (int)it.t.size()) {try {size_t n; v = std::stod(it.t[ig+3], &n); if(n==it.t[ig+3].size()) sp.dt.T = v;} catch(...) {}}
                }
                else it.fail("generator must be none, constant <T>, slip <slope> <rpm> or k_omega2 <K> [<T_max>]");
            }
            sp.pitch_joints = it.has("pitch_joints");
            sp.pitch_freq = it.get("pitch_freq", sp.pitch_freq);
            sp.pitch_zeta = it.get("pitch_zeta", sp.pitch_zeta);
            sp.blade_elements = (int)it.get("blade_elements", 0.0);
            sp.blade_beta = it.get("blade_beta", sp.blade_beta);
            sp.ea_per_mass = it.get("ea_per_mass", sp.ea_per_mass);
            sp.gj_factor = it.get("gj_factor", sp.gj_factor);
            if(sp.blade_elements < 0) it.fail("blade_elements must not be negative");
            sp.aero.rho = rho_air;
            sp.aero.skew = !it.has("no_skew");
            sp.aero.bem.tip_loss = !it.has("no_tip_loss");
            sp.aero.bem.hub_loss = !it.has("no_hub_loss");
            sp.aero.bem.tangential = !it.has("no_tangential");
            sp.aero.bem.induction = !it.has("no_induction");
            sp.aero.bem.drag_in_induction = !it.has("no_drag_induction");
            const int it_ = it.find("tower");
            if(it_ >= 0)
            {
                sp.aero.tower = true;
                sp.aero.tower_base = it.vec(it_+1);
                sp.aero.tower_top = it.vec(it_+4);
                sp.aero.D_base = it.num(it_+7);
                sp.aero.D_top = it.num(it_+8);
                if(it.has("tower_body")) sp.aero.tower_body = body_of(it, it.word("tower_body"));
            }
            try
            {
                const rotor_handles h = build_rotor(m, wind, sp);
                rotors[sp.name] = h;
                rotor_specs_[sp.name] = sp;
                joints[sp.name + "_shaft"] = h.shaft;
                if(h.gen_joint) joints[sp.name + "_gen_shaft"] = h.gen_joint;
                for(size_t q=0; q<h.pitch_joints.size(); ++q) joints[sp.name + "_pitch" + std::to_string(q + 1)] = h.pitch_joints[q];
            }
            catch(const std::exception &e) {it.fail(e.what());}
        }
        else if(k=="controller")
        {
            const std::string rn = it.word("rotor");
            if(!rotors.count(rn)) it.fail("unknown rotor '" + rn + "' (the rotor line must come first)");
            rotor_handles &h = rotors[rn];
            const rotor_spec &sp = rotor_specs_[rn];
            const double N = sp.dt.N, rpm = M_PI/30.0;
            turbine_controller_settings cs;
            if(it.has("tune"))
            {
                const double Pr = it.need("P_rated"), Or = it.need("rpm_rated")*rpm;
                rotor_def rd;
                rd.B = sp.B; rd.R_hub = sp.R_hub; rd.R_tip = sp.stations.back().r; rd.stations = sp.stations; rd.airfoils = sp.airfoils;
                const rotor_tuning tu = tune_rotor(rd, rho_air, Pr, Or, 25.0, sp.aero.bem);
                if(tu.K*Or*Or*Or >= Pr)
                    it.fail("the k omega^2 curve reaches the rated power below the rated speed (K Omega_rated^3 = "
                            + std::to_string(tu.K*Or*Or*Or) + " W at lambda_opt): lower rpm_rated or raise P_rated "
                            "(the NREL 5 MW: 0.82 of the rated power at the rated speed)");
                if(it.has("write") && !tu.write(path(it.word("write")))) it.fail("cannot write " + it.word("write"));
                const double wref = N*Or;
                cs.w_ref = wref;
                cs.w_rated = 0.99*wref;
                cs.w_cutin = 0.571*wref;
                cs.w_rgn2 = 0.742*wref;
                cs.P_rated = Pr;
                cs.K = tu.K/(N*N*N);
                cs.T_max = 1.1*Pr/cs.w_rated;
                cs.T_rate = 0.35*Pr/cs.w_rated;
                cs.theta_K = tu.theta_K;
                const double J = h.J + sp.J_gen*N*N;
                turbine_controller_settings::pi_gains(J, Or, N, tu.dPdtheta0, it.get("omega_n", 0.6), it.get("zeta", 0.7), cs.KP, cs.KI);
            }
            cs.K = it.get("K", cs.K);
            cs.w_cutin = it.get("rpm_cutin", cs.w_cutin/rpm)*rpm;
            cs.w_rgn2 = it.get("rpm_rgn2", cs.w_rgn2/rpm)*rpm;
            cs.w_rated = it.get("rpm_rated_torque", cs.w_rated/rpm)*rpm;
            if(!it.has("tune")) cs.P_rated = it.get("P_rated", cs.P_rated);
            cs.slip = it.get("slip", cs.slip);
            cs.T_max = it.get("T_max", cs.T_max);
            cs.T_rate = it.get("T_rate", cs.T_rate);
            cs.pitch_rgn3 = it.get("pitch_rgn3", cs.pitch_rgn3/deg)*deg;
            cs.w_ref = it.get("rpm_ref", cs.w_ref/rpm)*rpm;
            cs.KP = it.get("KP", cs.KP);
            cs.KI = it.get("KI", cs.KI);
            cs.theta_K = it.get("theta_K", cs.theta_K/deg)*deg;
            cs.theta_min = it.get("theta_min", cs.theta_min/deg)*deg;
            cs.theta_max = it.get("theta_max", cs.theta_max/deg)*deg;
            cs.theta_rate = it.get("theta_rate", cs.theta_rate/deg)*deg;
            cs.omega_c = it.get("corner", cs.omega_c/(2.0*M_PI))*2.0*M_PI;
            cs.dt = it.get("dt", 0.0);
            cs.pitch_control = !it.has("no_pitch");
            cs.torque_control = !it.has("no_torque");
            controller_load::io f;
            drivetrain *dt = h.drive;
            if(cs.torque_control) dt->settings.law = generator_law::command;
            f.speed = [dt](const model &mm) {return dt->omega_hss(mm);};
            f.torque = [dt](const model &mm) {return dt->generator_torque(mm);};
            f.set_torque = [dt](double T) {dt->set_torque_command(T);};
            if(h.pitch)
            {
                pitch_actuator *pa = h.pitch;
                f.pitch = [pa](const model &mm) {return pa->angle(mm, 0);};
                f.set_pitch = [pa](double th) {pa->set_command_all(th);};
            }
            else
            {
                rotor_aero *ra = h.aero;
                const int B = sp.B;
                f.pitch = [ra](const model&) {return ra->pitch(0, 0.0);};
                f.set_pitch = [ra, B](double th) {ra->set_pitch_table({}, {}); for(int q=0; q<B; ++q) ra->set_pitch(q, th);};
            }
            m.add_load(std::make_unique<controller_load>(cs, f, it.word(1)));
        }
        else if(k=="wind_probe")
        {
            if(!wind) it.fail("wind_probe needs a wind inflow (wind ...)");
            const int b = it.has("body") ? body_of(it, it.word("body")) : -1;
            m.add_load(std::make_unique<wind_probe>(m, wind, b, it.vec(2), it.word(1)));
        }
        else if(k=="fluid_coupling")
        {
            fluid_order = (int)it.get("order", (double)fluid_order);
            if(fluid_order < 0 || fluid_order > 2) it.fail("order must be 0, 1 or 2");
            if(it.has("no_predictor")) fluid_predictor = false;
            if(it.has("origin")) fluid_origin = it.need_vec("origin");
            if(it.has("record")) fluid_record = it.word("record");
        }
        else if(k!="gravity" && k!="rho" && k!="time" && k!="output" && k!="waves" && k!="body" && k!="joint"
                && k!="section" && k!="member" && k!="point_mass" && k!="static" && k!="modes"
                && k!="superelement" && k!="soil" && k!="soil_layer" && k!="pile" && k!="force" && k!="fluid_kinematics"
                && k!="point" && k!="line_type" && k!="bathymetry" && k!="wind" && k!="rho_air")
            it.fail("unknown keyword");
    }
    if(!fluid_bodies.empty())
    {
        std::vector<int> b;
        for(const auto &o : fluid_bodies)
            b.push_back(o.body);
        fluid = &static_cast<fluid_load&>(m.add_load(std::make_unique<fluid_load>(m, b, "fluid")));
        fluid->order = fluid_order;
        for(int k=0; k<(int)fluid_bodies.size(); ++k)
            if(fluid_bodies[k].added_mass.any())
                fluid->set_artificial_added_mass(k, fluid_bodies[k].added_mass.asDiagonal());
    }

    // initial displacements (the joints, the reference positions of the loads and the line attachments
    // belong to the configuration above), then the dynamic lines settle in the displaced configuration
    for(const auto &o : offsets)
        m.bodies[o.first].x += o.second;
    for(const auto &st : settles)
        st.first->settle(m, st.second);
}

}
