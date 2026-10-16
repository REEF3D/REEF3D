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
#include"coral/core/coral_rotation.h"
#include"coral/structure/coral_spring.h"
#include"coral/structure/coral_linear6.h"
#include"coral/hydro/coral_hydro_db.h"
#include"coral/hydro/coral_hydro_load.h"
#include"coral/hydro/coral_morison.h"
#include"coral/pto/coral_pto.h"
#include"coral/mooring/coral_mooring.h"
#include<cmath>
#include<fstream>
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

        [[noreturn]] void fail(const std::string &msg) const
        {
            throw std::runtime_error("CORAL input, line " + std::to_string(line) + " (" + t[0] + "): " + msg);
        }
        int find(const std::string &key, int from=1) const
        {
            for(int i=from; i<(int)t.size(); ++i)
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
    for(const auto &l : lines)
        items.push_back(item{l.first, l.second});
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
    std::vector<std::pair<mooring_line*, double>> settles;
    bool gravity_on = true;
    double wdepth = -1.0, wramp = 0.0;
    std::vector<item> wave_items;
    for(const item &it : items)
    {
        const std::string &k = it.t[0];
        if(k=="gravity")
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
            run.out_dt = it.get("dt", 0.0);
            if(it.has("file")) run.out_file = it.word("file");
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
            if(il >= 0) {s.control = pto_control::latching; s.t_hold = it.num(il+1); s.c_latch = it.num(il+2);}
            const int idc = it.find("declutching");
            if(idc >= 0) {s.control = pto_control::declutching; s.t_hold = it.num(idc+1);}
            m.add_load(std::make_unique<pto_load>(*joints[jn], it.has("rotational"), s, name));
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
            const int ik = it.find("k"), ic = it.find("c");
            for(int d=0; d<6; ++d)
            {
                if(ik >= 0) K(d,d) = it.num(ik+1+d);
                if(ic >= 0) C(d,d) = it.num(ic+1+d);
            }
            m.add_load(std::make_unique<linear6>(m, b, it.need_vec("point"), K, C, Vec6::Zero(), it.word(1)));
        }
        else if(k=="mooring")
        {
            const std::string name = it.word(1), type = it.word(2);
            const int ia = it.find("a"), ib = it.find("b");
            if(ia < 0 || ib < 0) it.fail("'a <body> <x y z>' and 'b <body> <x y z>' are required");
            mooring_end ea, eb;
            ea.body = body_of(it, it.word(ia+1)); ea.p = it.vec(ia+2);
            eb.body = body_of(it, it.word(ib+1)); eb.p = it.vec(ib+2);
            mooring_props p;
            p.L = it.need("L"); p.d = it.need("d"); p.mass = it.need("mass"); p.EA = it.need("EA");
            p.BA = it.get("BA", 0.0); p.CB = it.get("CB", 0.0);
            p.Cdn = it.get("Cdn", p.Cdn); p.Cdt = it.get("Cdt", p.Cdt); p.Can = it.get("Can", p.Can); p.Cat = it.get("Cat", p.Cat);
            p.kb = it.get("kb", p.kb); p.cb = it.get("cb", p.cb); p.mu = it.get("mu", 0.0);
            p.rho = rho;
            if(type=="quasi_static")
                m.add_load(std::make_unique<mooring_qs>(m, ea, eb, p, name));
            else if(type=="dynamic")
            {
                const double zbed = it.get("seabed", std::min(ea.p(2), eb.p(2)));
                auto l = std::make_unique<mooring_line>(m, ea, eb, p, (int)it.get("nseg", 20), waves, it.vec("current", Vec3::Zero()),
                                                        zbed, it.get("dt_line", 0.0), name);
                if(it.has("settle"))
                    settles.push_back({l.get(), it.need("settle")});
                m.add_load(std::move(l));
            }
            else
                it.fail("mooring type must be quasi_static or dynamic");
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
            m.add_load(std::make_unique<morison_load>(m, b, std::vector<morison_member>{mb}, waves, rho, it.vec("current", Vec3::Zero()), it.word(1)));
        }
        else if(k!="gravity" && k!="rho" && k!="time" && k!="output" && k!="waves" && k!="body" && k!="joint")
            it.fail("unknown keyword");
    }

    // initial displacements (the joints, the reference positions of the loads and the line attachments
    // belong to the configuration above), then the dynamic lines settle in the displaced configuration
    for(const auto &o : offsets)
        m.bodies[o.first].x += o.second;
    for(const auto &st : settles)
        st.first->settle(m, st.second);
}

}
