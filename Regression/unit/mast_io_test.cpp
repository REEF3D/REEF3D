// Standalone tests of the REEF3D::MAST input and output (src/mast/io): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         mast_io_test.cpp ../../src/mast/core/*.cpp ../../src/mast/structure/*.cpp ../../src/mast/hydro/*.cpp
//         ../../src/mast/pto/*.cpp ../../src/mast/mooring/*.cpp ../../src/mast/io/*.cpp -o mast_io_test
// Run:    ./mast_io_test
//
// Tests: a model with every keyword (bodies with inertia products, Euler angles and velocities,
// joints, waves, hydrodynamic database, PTO, spring, linear6, quasi-static and dynamic mooring,
// Morison) is read with the right numbers; line continuation and comments; errors name the line;
// a pendulum read from the input runs bit for bit like the same model built in code; the output
// file has one column per channel.
// Architect: Hans Bihs
#include"mast/core/mast_gravity.h"
#include"mast/core/mast_constraint.h"
#include"mast/core/mast_timeint.h"
#include"mast/hydro/mast_hydro_db.h"
#include"mast/io/mast_input.h"
#include"mast/io/mast_output.h"
#include<cmath>
#include<cstdio>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>

using namespace mast;

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
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
    db.write("/tmp/mast_io_test.mhd");

    const std::string text = R"(
# test model
gravity 0 0 -9.81
rho 1000
time dt 0.01 end 1.0 rho_inf 0.7
output dt 0.05 file /tmp/mast_io_test_out.dat
body buoy mass 9.0 inertia 0.1 0.2 0.3 0.01 0.02 0.03 cog 0 0 -0.05 \
     rpy 0 0 90 velocity 0 0 0.1
body plate mass 2.0 inertia 1 1 1 cog 0 0 -0.5
joint heave prismatic ground buoy point 0 0 0 axis 0 0 1
joint hinge revolute buoy plate point 0 0 -0.3 axis 0 1 0
waves depth 1.0
waves ramp 2.0
waves regular H 0.04 T 1.0
waves jonswap Hs 0.05 Tp 1.2 gamma 3.3 wmin 3 wmax 11 n 16 seed 4
hydro bem file /tmp/mast_io_test.mhd bodies buoy damping 3 40 drag 3 5
pto p1 joint heave c 40 stops -0.1 0.1 1e5
pto p2 joint hinge rotational k 3 c 0.5
spring tether a ground 0 0 -1 b buoy 0 0 -0.15 k 150 T0 15 tension_only
linear6 moor body plate point 0 0 -0.5 k 1 2 3 4 5 6 c 0 0 0 0 0 0
mooring m1 quasi_static a ground 10 0 -1 b buoy 0 0 -0.1 L 11 d 0.01 mass 0.5 EA 1e6
mooring m2 dynamic a ground -10 0 -1 b buoy 0 0 -0.1 L 11 d 0.01 mass 0.5 EA 1e6 nseg 8
morison leg body plate p1 0 0 -0.8 p2 0 0 -0.4 D 0.1 Cd 1.0 Ca 1.0 n 4
)";
    mast_case c;
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
    check(c.run.dt==0.01 && c.run.t_end==1.0 && c.run.ti.rho_inf==0.7 && c.run.out_dt==0.05 && c.run.out_file=="/tmp/mast_io_test_out.dat", "time and output settings");
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
        try {mast_case c; c.read_string(text);} catch(const std::exception &e) {return std::string(e.what());}
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
    mast_case c;
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

int main()
{
    test_full_model();
    test_errors();
    test_same_as_code();
    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
