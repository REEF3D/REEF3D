// Standalone tests of the REEF3D::CORAL input and output (src/coral/io): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         coral_io_test.cpp ../../src/coral/core/*.cpp ../../src/coral/structure/*.cpp ../../src/coral/hydro/*.cpp
//         ../../src/coral/pto/*.cpp ../../src/coral/mooring/*.cpp ../../src/coral/io/*.cpp -o coral_io_test
// Run:    ./coral_io_test
//
// Tests: a model with every keyword (bodies with inertia products, Euler angles and velocities,
// joints, waves, hydrodynamic database, PTO, spring, linear6, quasi-static and dynamic mooring,
// Morison) is read with the right numbers; line continuation and comments; errors name the line;
// a pendulum read from the input runs bit for bit like the same model built in code; the output
// file has one column per channel; structure keywords (section, member, point_mass, static, modes,
// output vtp) build the same frame as the code and give bit for bit the same modes, the mode table
// (kinetic energy split) and the VTP file are written; offshore keywords (superelement, soil,
// soil_layer, pile, morison on frame and superelement members, point_mass at an offset, force, linear6
// with full matrices) build the expected loads and the model finds its static equilibrium.
// Architect: Hans Bihs
#include"coral/core/coral_gravity.h"
#include"coral/core/coral_constraint.h"
#include"coral/core/coral_timeint.h"
#include"coral/hydro/coral_hydro_db.h"
#include"coral/io/coral_input.h"
#include"coral/io/coral_output.h"
#include"coral/io/coral_modes_out.h"
#include"coral/io/coral_vtp.h"
#include"coral/core/coral_static.h"
#include"coral/structure/coral_frame.h"
#include"coral/linear/coral_modes.h"
#include"coral/structure/coral_superelement.h"
#include"coral/structure/coral_linear6.h"
#include"coral/core/coral_force.h"
#include"coral/soil/coral_soil.h"
#include"coral/hydro/coral_member_morison.h"
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iostream>
#include<memory>
#include<sstream>
#include<string>

using namespace coral;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static std::string fmt_n(const char *f, int n)
{
    char b[256];
    std::snprintf(b, sizeof(b), f, n);
    return b;
}

static std::string fmt_modes(const mode_result &r)
{
    char b[256];
    std::snprintf(b, sizeof(b), "input file and code give bit for bit the same modes (%.4f %.4f %.4f %.4f Hz)",
                  r.frequency(0), r.frequency(1), r.frequency(2), r.frequency(3));
    return b;
}

static void test_full_model()
{
    std::cout<<"model with every keyword"<<std::endl;
    // a one-body heave database
    hydro_db db;
    db.resize(1, 3, 1);
    db.rho = 1000.0;
    db.w = {1.0, 2.0, 3.0};
    db.C(2,2) = 700.0;
    for(int k=0; k<3; ++k) {db.A[k](2,2) = 5.0; db.B[k](2,2) = 2.0*k;}
    db.body[0].volume = 0.01;
    db.write("/tmp/coral_io_test.chd");

    const std::string text = R"(
# test model
gravity 0 0 -9.81
rho 1000
time dt 0.01 end 1.0 rho_inf 0.7
output dt 0.05 file /tmp/coral_io_test_out.dat
body buoy mass 9.0 inertia 0.1 0.2 0.3 0.01 0.02 0.03 cog 0 0 -0.05 \
     rpy 0 0 90 velocity 0 0 0.1
body plate mass 2.0 inertia 1 1 1 cog 0 0 -0.5
joint heave prismatic ground buoy point 0 0 0 axis 0 0 1
joint hinge revolute buoy plate point 0 0 -0.3 axis 0 1 0
waves depth 1.0
waves ramp 2.0
waves regular H 0.04 T 1.0
waves jonswap Hs 0.05 Tp 1.2 gamma 3.3 wmin 3 wmax 11 n 16 seed 4
hydro bem file /tmp/coral_io_test.chd bodies buoy damping 3 40 drag 3 5
pto p1 joint heave c 40 stops -0.1 0.1 1e5
pto p2 joint hinge rotational k 3 c 0.5
spring tether a ground 0 0 -1 b buoy 0 0 -0.15 k 150 T0 15 tension_only
linear6 moor body plate point 0 0 -0.5 k 1 2 3 4 5 6 c 0 0 0 0 0 0
mooring m1 quasi_static a ground 10 0 -1 b buoy 0 0 -0.1 L 11 d 0.01 mass 0.5 EA 1e6
mooring m2 dynamic a ground -10 0 -1 b buoy 0 0 -0.1 L 11 d 0.01 mass 0.5 EA 1e6 nseg 8
morison leg body plate p1 0 0 -0.8 p2 0 0 -0.4 D 0.1 Cd 1.0 Ca 1.0 n 4
)";
    coral_case c;
    bool ok = true;
    try {c.read_string(text, "/tmp");} catch(const std::exception &e) {std::cout<<"  "<<e.what()<<std::endl; ok = false;}
    check(ok, "reads without error");
    if(!ok) return;
    check(c.m.nbody()==2 && c.m.bodies[0].name=="buoy" && c.m.bodies[1].name=="plate", "two bodies by name");
    const rigid_body &b = c.m.bodies[0];
    check(b.mass==9.0 && b.inertia(0,0)==0.1 && b.inertia(1,1)==0.2 && b.inertia(0,1)==0.01 && b.inertia(1,2)==0.03, "mass and inertia with products");
    check(std::fabs(b.R(0,1) + 1.0) < 1.0e-12 && std::fabs(b.v(2) - 0.1) < 1.0e-15 && std::fabs(b.x(2) + 0.05) < 1.0e-15, "Euler angles (yaw 90 deg), velocity and CoG (continued line)");
    check(c.joints.size()==2 && c.m.ncon()==10, "two joints, 10 constraint equations");
    check(c.waves && c.waves->components().size()==17 && c.waves->depth()==1.0 && std::fabs(c.waves->ramp(1.0) - 0.5) < 1.0e-15, "waves: 17 components, depth and ramp");
    // gravity, hydro, 2 pto, spring, linear6, 2 mooring, morison
    check(c.m.load_list.size()==9, "nine loads (gravity included)");
    check(c.run.dt==0.01 && c.run.t_end==1.0 && c.run.ti.rho_inf==0.7 && c.run.out_dt==0.05 && c.run.out_file=="/tmp/coral_io_test_out.dat", "time and output settings");
    // run it briefly
    timeint ti(c.run.ti);
    ti.initialize(c.m, 0.0);
    time_series out;
    out.open(c.run.out_file, c.m);
    for(int i=0; i<20; ++i) {ti.step(c.m, i*0.01, 0.01); out.write(c.m, (i+1)*0.01);}
    out.close();
    std::ifstream in(c.run.out_file);
    std::string l1, l2, l3;
    std::getline(in, l1); std::getline(in, l2); std::getline(in, l3);
    std::istringstream h(l2), d(l3);
    std::string w;
    int nh = -1, nd = 0;     // '#' is not a column
    while(h >> w) ++nh;
    while(d >> w) ++nd;
    check(nh==(int)out.columns().size() && nd==nh && ti.stats.failures==0, "output: header and values have one column per channel (" + std::to_string(nd) + ")");
}

static void test_errors()
{
    std::cout<<"errors name the line"<<std::endl;
    auto msg = [](const std::string &text)
    {
        try {coral_case c; c.read_string(text);} catch(const std::exception &e) {return std::string(e.what());}
        return std::string();
    };
    std::string e = msg("body a mass 1\n\nfoo 3\n");
    check(e.find("line 3")!=std::string::npos && e.find("unknown keyword")!=std::string::npos, "unknown keyword: " + e);
    e = msg("body a mass 1\njoint j revolute ground b point 0 0 0\n");
    check(e.find("line 2")!=std::string::npos && e.find("unknown body 'b'")!=std::string::npos, "unknown body: " + e);
    e = msg("body a inertia 1 1 1\n");
    check(e.find("line 1")!=std::string::npos && e.find("'mass' is required")!=std::string::npos, "missing value: " + e);
    e = msg("body a mass x1\n");
    check(e.find("is not a number")!=std::string::npos, "not a number: " + e);
}

static void test_same_as_code()
{
    std::cout<<"input file and code give the same run"<<std::endl;
    coral_case c;
    c.read_string("body bob mass 2 inertia 0.008 0.008 0.008 cog 1.5 0 0\njoint pin revolute ground bob point 0 0 0 axis 0 1 0\n");
    model m;
    rigid_body b;
    b.name = "bob";
    b.mass = 2.0;
    b.inertia = Vec3(0.008, 0.008, 0.008).asDiagonal();
    b.x = Vec3(1.5, 0, 0);
    m.add_body(b);
    m.add_load(std::make_unique<gravity_load>());
    m.add_constraint(std::make_unique<joint>(m, joint_type::revolute, -1, 0, Vec3::Zero(), Vec3::UnitY()));
    timeint t1, t2;
    t1.initialize(c.m, 0.0);
    t2.initialize(m, 0.0);
    for(int i=0; i<200; ++i) {t1.step(c.m, i*0.005, 0.005); t2.step(m, i*0.005, 0.005);}
    check(c.m.bodies[0].x==m.bodies[0].x && c.m.bodies[0].R==m.bodies[0].R, "pendulum: bit for bit the same state after 1 s");
}

static void test_structure()
{
    std::cout<<"structure keywords"<<std::endl;
    const std::string text = R"(
section a tube D 0.5 t 0.02 E 2.1e11 G 8.08e10 rho 7850
section b EA 1e9 EIy 2e7 EIz 3e7 GJ 1e7 m 100 rIx 1 GAy 4e8
member post p1 0 0 0 p2 0 0 10 n 4 section a orient 1 0 0 end1 ground
member arm p1 0 0 10 p2 5 0 10 n 2 section b orient 0 0 1
point_mass arm:2 mass 500 inertia 1 2 3
static
modes 4
output vtp 0.5
time dt 0.01 end 0
)";
    coral_case c;
    bool ok = true;
    try {c.read_string(text);} catch(const std::exception &e) {std::cout<<"  "<<e.what()<<std::endl; ok = false;}
    check(ok, "reads without error");
    if(!ok) return;
    const run_settings &r = c.run;
    check(r.static_ini && r.modes==4 && !r.modes_damped && r.vtp_dt==0.5 && r.t_end==0.0, "static, modes, output vtp, time");
    int tip = -1;
    for(int b=0; b<c.m.nbody(); ++b) if(c.m.bodies[b].name=="arm:2") tip = b;
    check(c.m.nbody()==6 && tip >= 0 && c.m.bodies[tip].mass > 500.0 && c.m.bodies[tip].inertia(2,2) > 3.0,
          "six node bodies (shared corner node), point mass and inertia on 'arm:2'");
    // the same frame in code
    model m;
    {
        frame_builder fb(m);
        member_spec post, arm;
        post.name = "post"; post.p2 = Vec3(0, 0, 10); post.n = 4; post.s1 = beam_section::tube(0.5, 0.02, 2.1e11, 8.08e10, 7850.0);
        post.orient = Vec3::UnitX(); post.body1 = -1;
        beam_section b;
        b.EA = 1e9; b.EIy = 2e7; b.EIz = 3e7; b.GJ = 1e7; b.m = 100.0; b.rIx = 1.0; b.GAy = 4e8;
        arm.name = "arm"; arm.p1 = Vec3(0, 0, 10); arm.p2 = Vec3(5, 0, 10); arm.n = 2; arm.s1 = b; arm.orient = Vec3::UnitZ();
        fb.add_member(post);
        fb.add_member(arm);
        fb.finish();
    }
    m.bodies[tip].mass += 500.0;
    m.bodies[tip].inertia += Mat3(Vec3(1, 2, 3).asDiagonal());
    m.add_load(std::make_unique<gravity_load>());     // after the members, as the input does
    const static_result s1 = static_solve(c.m, 0.0), s2 = static_solve(m, 0.0);
    check(s1.converged && s2.converged, "static equilibrium of both");
    const mode_result r1 = modes(c.m, 0.0, 4), r2 = modes(m, 0.0, 4);
    check(r1.ok && r2.ok && r1.frequency==r2.frequency, fmt_modes(r1));
    write_modes("/tmp/coral_io_test_", c.m, r1);
    std::ifstream mf("/tmp/coral_io_test_modes.dat");
    std::string line;
    int nl = 0;
    bool sum1 = true;
    while(std::getline(mf, line))
    {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream l(line);
        double q, f, T, z, e, es = 0.0;
        l>>q>>f>>T>>z;
        for(int k=0; k<6; ++k) {l>>e; es += e;}
        sum1 = sum1 && std::fabs(es - 1.0) < 1.0e-3 && std::fabs(f - r1.frequency(nl)) < 1.0e-6*f;
        ++nl;
    }
    check(nl==4 && sum1, "modes.dat: four modes, kinetic energy split sums to 1");
    std::ifstream sf("/tmp/coral_io_test_mode_shapes.dat");
    int ns = 0;
    while(std::getline(sf, line)) if(!line.empty() && line[0]!='#') ++ns;
    check(ns==4*6, "mode_shapes.dat: one line per body and mode");
    check(write_vtp(c.m, "/tmp/coral_io_test.vtp"), "VTP written");
    std::ifstream vf("/tmp/coral_io_test.vtp");
    std::stringstream vs;
    vs<<vf.rdbuf();
    check(vs.str().find("NumberOfLines=\"6\"")!=std::string::npos, "VTP: six beam elements as lines");
}

static void test_offshore()
{
    std::cout<<"offshore keywords"<<std::endl;
    const std::string text = R"(
rho 1025
waves depth 30
waves regular H 2 T 9
section leg tube D 1.2 t 0.04 E 2.1e11 G 8.08e10 rho 7850
section tw tube D 5 t 0.04 E 2.1e11 G 8.08e10 rho 7850
section pl tube D 2 t 0.05 E 2.1e11 G 8.08e10 rho 7850
body tp mass 2e5 inertia 4e6 4e6 6e6 cog 0 0 10
member l1 p1 8 0 -30 p2 3 0 10 n 6 section leg end1 ground end2 tp orient 0 0 1
member l2 p1 -4 6.93 -30 p2 -1.5 2.6 10 n 6 section leg end1 ground end2 tp orient 0 0 1
member l3 p1 -4 -6.93 -30 p2 -1.5 -2.6 10 n 6 section leg end1 ground end2 tp orient 0 0 1
superelement jacket members l1 l2 l3 modes 12 zeta 0.01
member tower p1 0 0 10 p2 0 0 70 n 10 section tw end1 tp orient 1 0 0
point_mass tp mass 1e4 at 3 0 10
member p p1 20 0 -30 p2 20 0 -60 n 15 section pl orient 1 0 0
soil s1 mudline -30
soil_layer s1 clay top 0 bottom 20 su 20e3 su_bottom 80e3 gamma 8e3 eps50 0.01
soil_layer s1 sand top 20 bottom 40 phi 33 gamma 10e3 k 2e7
pile p soil s1 D 2
morison wl1 member l1 Cd 1.0 Ca 1.0
morison wt member tower Cd 0.7 Ca 1.0
force push body tp point 0 0 12 F 1e5 0 0 M 0 0 1e4 ramp 2
linear6 cs body tp point 0 0 10 K 1e9 0 0 0 0 0 0 1e9 0 0 0 0 0 0 1e9 0 0 0 0 0 0 1e10 0 0 0 0 0 0 1e10 0 0 0 0 0 0 1e10
static
time dt 0.05 end 0
)";
    coral_case c;
    bool ok = true;
    try {c.read_string(text);} catch(const std::exception &e) {std::cout<<"  "<<e.what()<<std::endl; ok = false;}
    check(ok, "reads without error");
    if(!ok) return;
    int nse = 0, nsoil = 0, nmor = 0, nl6 = 0, nmodes = 0, nforce = 0;
    for(const auto &l : c.m.load_list)
    {
        if(const auto *se = dynamic_cast<const superelement*>(l.get())) {++nse; nmodes = se->nmode();}
        if(dynamic_cast<const soil_spring*>(l.get())) ++nsoil;
        if(dynamic_cast<const member_morison*>(l.get())) ++nmor;
        if(dynamic_cast<const linear6*>(l.get())) ++nl6;
        if(dynamic_cast<const point_force*>(l.get())) ++nforce;
    }
    check(nse==1 && nmodes==12 && c.m.ncoord()==12, "superelement with 12 modes (12 generalized coordinates)");
    check(nforce==1, "point force");
    check(nsoil==15 && nmor==2 && nl6==1 && c.m.point_masses.size()==2, fmt_n("15 pile springs below the mudline (%d), two Morison members, linear6, point masses at offsets (point_mass, tower end on tp)", nsoil));
    // tower nodes and the pile are frame nodes; the jacket legs are not
    check(c.m.find_body("tower:3") >= 0 && c.m.find_body("p:3") >= 0 && c.m.find_body("l1:3") < 0, "frame nodes for the tower and the pile, none for the superelement legs");
    const static_result r = static_solve(c.m, 0.0);
    check(r.converged, fmt_n("static equilibrium in %d iterations", r.iterations));
    check(write_vtp(c.m, "/tmp/coral_io_test_offshore.vtp"), "VTP written");
    std::ifstream vf("/tmp/coral_io_test_offshore.vtp");
    std::stringstream vs;
    vs<<vf.rdbuf();
    check(vs.str().find("\n3\n")!=std::string::npos, "VTP: superelement elements (kind 3)");
}

int main()
{
    test_full_model();
    test_errors();
    test_same_as_code();
    test_structure();
    test_offshore();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
