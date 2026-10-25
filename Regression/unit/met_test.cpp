// Architect: Hans Bihs
// Standalone verification of the REEF3D atmospheric forcing kernels (W 200): date-times, drag laws (met_drag.h),
// geographic frame (met_frame) against tools/seastate_forcing.py, binary met file reader (met_file): header,
// file written by tools/met_forcing.py (met_test.r3m), bilinear/linear interpolation exact for linear fields,
// windows (only the needed part read, the same values), record advance, clamping in space and time, uniform
// series, moving grids, bad files; uniform and file sources (met_source.h). No MPI, no REEF3D binary.
// Build:  g++ -O2 -std=c++20 -I../../src met_test.cpp ../../src/met_file.cpp ../../src/met_frame.cpp -o met_test
// Run:    ./met_test
#include"met_drag.h"
#include"met_file.h"
#include"met_frame.h"
#include"met_source.h"
#include<cmath>
#include<cstdint>
#include<cstdio>
#include<cstring>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>

static int nfail = 0;

static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static const double pi = 3.14159265358979323846;

// ---------------------------------------------------------------------------------------------
// writer of the format (as tools/met_forcing.py write_r3m)
// ---------------------------------------------------------------------------------------------
struct field
{
    // u, v, pa linear in x, y and t (t: s after the first record), the coefficients of met_forcing.py ANALYTIC
    double c[12] = {10.0, 1.0e-4, -2.0e-4, 1.0e-4, -3.0, 2.0e-4, 1.0e-4, -5.0e-5, 101000.0, -0.01, 0.02, -0.005};

    void at(double x, double y, double t, double &u, double &v, double &pa) const
    {
        u  = c[0] + c[1]*x + c[2]*y  + c[3]*t;
        v  = c[4] + c[5]*x + c[6]*y  + c[7]*t;
        pa = c[8] + c[9]*x + c[10]*y + c[11]*t;
    }
};

static void write_file(const std::string &name, int nx, int ny, double x0, double y0, double dx, double dy,
                       const std::vector<double> &times, const field &F, int nvar=3, bool moving=false,
                       const std::vector<double> &ox = {}, const std::vector<double> &oy = {})
{
    std::ofstream out(name.c_str(), std::ios::binary);
    out.write("REEF3DMT",8);
    int32_t iv[8] = {1, nx, ny, nvar, moving ? 1 : 0, int32_t(times.size()), 1, 0};
    out.write((char*)iv,sizeof(iv));
    double dv[13] = {x0, y0, dx, dy, 5.9, 62.3, 20.0, 0.0, 0.0, 600.0, 0.0, 0.0, 0.0};
    out.write((char*)dv,sizeof(dv));

    for(size_t n=0; n<times.size(); ++n)
    {
        const double xr = moving ? ox[n] : x0, yr = moving ? oy[n] : y0;
        double h[4] = {times[n], xr, yr, 0.0};
        out.write((char*)h,sizeof(h));

        std::vector<float> a[3];
        for(int q=0; q<3; ++q)
        a[q].resize(size_t(nx)*ny);

        for(int j=0; j<ny; ++j)
        for(int i=0; i<nx; ++i)
        {
            double u, v, pa;
            F.at(xr + i*dx, yr + j*dy, times[n]-times[0], u, v, pa);
            a[0][j*nx+i] = float(u);
            a[1][j*nx+i] = float(v);
            a[2][j*nx+i] = float(pa);
        }

        for(int q=0; q<nvar; ++q)
        out.write((char*)a[q].data(), 4*a[q].size());
    }
}

// ---------------------------------------------------------------------------------------------
static void test_datetime()
{
    std::cout<<"date-times"<<std::endl;
    check(met_datetime(19700101.0)==0.0, "19700101.000000 -> 0");
    check(met_datetime(20181010.173000)==1539192600.0, "20181010.173000 -> 1539192600 (Python datetime)");
    check(met_datetime(20171231.235959)==1514764799.0, "20171231.235959 -> 1514764799");
    check(met_datetime(20240229.120000)==1709208000.0, "20240229.120000 (leap day) -> 1709208000");
    check(met_date(1539192600.0)=="2018-10-10 17:30:00", "met_date back to 2018-10-10 17:30:00");
}

static void test_drag()
{
    std::cout<<"drag laws (as wind_f, A 570)"<<std::endl;
    bool ok = true;

    for(double U : {2.0, 7.5, 10.0, 20.0, 35.0, 41.0, 60.0})
    {
        double g = U<41.0 ? 0.001 * (0.75+0.067*U) : 0.001 * 3.5;
        double w = U<7.5 ? 1.2875e-3 : 0.001 * (0.8 + 0.065*U);
        double s = 0.001 * (0.63+0.066*U);
        double z = 0.001 * (0.55 + 2.97*(U/31.5) - 1.49*std::pow(U/31.5,2.0));
        ok = ok && met_drag(1,U)==g && met_drag(2,U)==w && met_drag(3,U)==s && std::fabs(met_drag(4,U)-std::max(z,0.0))<1.0e-18;
    }
    check(ok, "Garratt, Wu, Smith and Banke bitwise, Zijlema to 1e-18 at 2 - 60 m/s");
    check(std::fabs(met_drag(1,20.0)-2.09e-3)<1.0e-15, "Garratt Cd(20 m/s) = 2.09e-3");
    check(met_drag(4,80.0)==0.0, "Zijlema not negative above 68 m/s");
    check(std::fabs(met_drag(4,31.5)-2.03e-3)<1.0e-15, "Zijlema Cd(Uref) = 2.03e-3");
}

static void test_frame()
{
    std::cout<<"frame (transverse Mercator, tools/seastate_forcing.py --origin --rotation)"<<std::endl;

    // lon0 lat0 rot lon lat -> x y north, from seastate_forcing.Frame
    const double ref[6][8] = {
    {5.9,62.3,20.0,6.6000000000000005,61.9,19305.6393789981,-54148.9097183036,1.232507813161},
    {5.9,62.3,20.0,3.9000000000000004,63.8,-34681.3615721215,191756.0563017148,1.190407757973},
    {-85.5,29.5,0.0,-84.8,29.1,68012.2183366765,-44275.9099507434,1.576738269331},
    {-85.5,29.5,0.0,-87.5,31.0,-190643.4826203823,168506.5350760436,1.552812713695},
    {150.0,-20.0,-35.0,150.7,-20.4,85362.7398963341,5284.4089466093,2.177402777839},
    {150.0,-20.0,-35.0,148.0,-18.5,-267783.1653065431,14685.6763005126,2.192741605452}};

    double ex = 0.0, en = 0.0, eg = 0.0;

    for(auto &r : ref)
    {
        met_frame f(r[0],r[1],r[2],0.0,0.0);
        double x, y, lon, lat;
        f.to_model(r[3],r[4],x,y);
        ex = std::max(ex, std::max(std::fabs(x-r[5]), std::fabs(y-r[6])));
        en = std::max(en, std::fabs(f.north(x,y)-r[7]));
        f.to_geo(x,y,lon,lat);
        eg = std::max(eg, std::max(std::fabs(lon-r[3]), std::fabs(lat-r[4])));
    }
    check(ex<1.0e-6, "model coordinates = Python to 1e-6 m (max "+std::to_string(ex)+")");
    check(en<1.0e-8, "north = Python to 1e-8 rad");
    check(eg<1.0e-11, "to_geo(to_model) = identity to 1e-11 deg");

    // offset: X0, Y0 is where (lon0, lat0) lands
    met_frame f(5.9,62.3,20.0,1000.0,-2000.0);
    double x, y;
    f.to_model(5.9,62.3,x,y);
    check(std::fabs(x-1000.0)<1.0e-9 && std::fabs(y+2000.0)<1.0e-9, "W 206: (lon0, lat0) at (X0, Y0)");
    met_frame g(10.0,0.0,0.0,0.0,0.0);
    check(std::fabs(g.north(0.0,0.0)-0.5*pi)<1.0e-9, "rot 0: north along +y");
    met_frame h(10.0,0.0,30.0,0.0,0.0);
    check(std::fabs(h.north(0.0,0.0)-(0.5*pi-pi/6.0))<1.0e-9, "rot 30: north at 60 deg from the model x axis");
}

static void test_python_file()
{
    std::cout<<"met file written by tools/met_forcing.py (met_test.r3m)"<<std::endl;
    met_file m;
    std::string err;
    const bool ok = m.open("met_test.r3m",err);
    check(ok, "open: "+(ok ? std::string("ok") : err));
    if(!ok) return;

    const double t0 = met_datetime(20181010.120000);
    check(m.nx==5 && m.ny==4 && m.nvar==3 && m.nrec==3 && m.moving==0, "header 5 x 4 nodes, 3 variables, 3 records");
    check(m.dx==1000.0 && m.dy==1000.0 && m.x0==0.0 && m.y0==0.0, "grid spacing and origin");
    check(m.time(0)==t0 && m.time(2)==t0+7200.0, "record times 2018-10-10 12:00 - 14:00");

    field F;
    double err_u = 0.0, err_p = 0.0;

    for(double t : {0.0, 1234.5, 3600.0, 5000.0, 7200.0})
    {
        m.advance(t0+t,err);
        for(double x : {0.0, 333.0, 2500.0, 3999.0, 4000.0})
        for(double y : {0.0, 1500.0, 2999.9})
        {
            double u, v, pa, ue, ve, pe;
            m.at(t0+t,x,y,u,v,pa);
            F.at(x,y,t,ue,ve,pe);
            err_u = std::max(err_u, std::max(std::fabs(u-ue), std::fabs(v-ve)));
            err_p = std::max(err_p, std::fabs(pa-pe));
        }
    }
    check(err_u<2.0e-6, "u10, v10 = the linear field (bilinear in space, linear in time; float32 nodes)");
    check(err_p<0.02,   "pa = the linear field to the float32 resolution");
}

static void test_window()
{
    std::cout<<"windows, advance, clamping"<<std::endl;
    field F;
    std::vector<double> times;
    for(int n=0; n<6; ++n)
    times.push_back(1.0e9 + 3600.0*n);

    write_file("met_test_tmp.r3m",101,81,-5000.0,2000.0,500.0,500.0,times,F);

    met_file a, b;
    std::string err;
    check(a.open("met_test_tmp.r3m",err) && b.open("met_test_tmp.r3m",err), "open the written file");

    // b reads the box 10 - 20 km x 10 - 15 km only
    b.window(10000.0,20000.0,10000.0,15000.0);

    bool same = true;
    for(double t : {0.0, 1800.0, 3600.0, 9000.0, 17999.0})
    {
        a.advance(1.0e9+t,err);
        b.advance(1.0e9+t,err);
        for(double x=10000.0; x<=20000.0; x+=777.0)
        for(double y=10000.0; y<=15000.0; y+=333.0)
        {
            double u1,v1,p1,u2,v2,p2;
            a.at(1.0e9+t,x,y,u1,v1,p1);
            b.at(1.0e9+t,x,y,u2,v2,p2);
            same = same && u1==u2 && v1==v2 && p1==p2;
        }
    }
    check(same, "the window gives bitwise the values of the full grid inside the box");
    check(b.values_read*10<a.values_read, "the window reads less than 10 % of the full grid ("+std::to_string(b.values_read)+" of "+std::to_string(a.values_read)+" values)");

    // advancing record by record reads one record per step
    met_file c;
    c.open("met_test_tmp.r3m",err);
    c.window(0.0,1000.0,3000.0,4000.0);
    c.advance(1.0e9,err);
    const long long r0 = c.values_read;
    c.advance(1.0e9+3600.0,err);
    const long long r1 = c.values_read;
    c.advance(1.0e9+3700.0,err);
    check(r1==r0+r0/2 && c.values_read==r1, "the next record only is read when the time passes a record, none within");

    // a grown window reads the records again, with the new box
    c.window(0.0,30000.0,3000.0,4000.0);
    c.advance(1.0e9+3700.0,err);
    double u,v,pa,ue,ve,pe;
    c.at(1.0e9+3700.0,25000.0,3500.0,u,v,pa);
    F.at(25000.0,3500.0,3700.0,ue,ve,pe);
    check(c.values_read>r1 && std::fabs(u-ue)<1.0e-5, "a grown window: records read again, values in the new part");

    // clamping in time and space
    a.advance(1.0e9-5000.0,err);
    a.at(1.0e9-5000.0,0.0,5000.0,u,v,pa);
    F.at(0.0,5000.0,0.0,ue,ve,pe);
    check(std::fabs(u-ue)<1.0e-5, "before the first record: the first record");
    a.advance(1.0e9+1.0e6,err);
    a.at(1.0e9+1.0e6,0.0,5000.0,u,v,pa);
    F.at(0.0,5000.0,18000.0,ue,ve,pe);
    check(std::fabs(u-ue)<1.0e-5, "after the last record: the last record");
    a.at(1.0e9+1.0e6,-1.0e6,5000.0,u,v,pa);
    F.at(-5000.0,5000.0,18000.0,ue,ve,pe);
    check(std::fabs(u-ue)<1.0e-5, "outside the grid: the value at the edge");
    std::remove("met_test_tmp.r3m");
}

static void test_uniform_and_moving()
{
    std::cout<<"uniform series, moving grid, sources"<<std::endl;
    field F;
    std::vector<double> times = {0.0, 600.0, 1200.0};
    write_file("met_test_tmp.r3m",1,1,0.0,0.0,1.0,1.0,times,F,2);

    met_file m;
    std::string err;
    check(m.open("met_test_tmp.r3m",err), "nx = ny = 1, two variables");
    m.advance(900.0,err);
    double u,v,pa,ue,ve,pe;
    m.at(900.0,12345.0,-777.0,u,v,pa);
    F.at(0.0,0.0,900.0,ue,ve,pe);
    check(std::fabs(u-ue)<1.0e-6 && std::fabs(v-ve)<1.0e-6 && pa==101325.0, "uniform in space, linear in time, pa 101325 Pa without the pressure");

    // moving: the origin of every record moves 2 km east per record; the field is fixed in space
    std::vector<double> ox = {0.0, 2000.0, 4000.0}, oy = {0.0, 500.0, 1000.0};
    write_file("met_test_tmp.r3m",21,21,0.0,0.0,250.0,250.0,times,F,3,true,ox,oy);
    met_file mv;
    check(mv.open("met_test_tmp.r3m",err) && mv.moving==1, "moving grid");
    mv.window(4500.0,5000.0,1500.0,2000.0);
    mv.advance(300.0,err);
    mv.at(300.0,4600.0,1700.0,u,v,pa);
    F.at(4600.0,1700.0,300.0,ue,ve,pe);
    check(std::fabs(u-ue)<1.0e-5 && std::fabs(pa-pe)<0.02, "moving grid: each record at its own origin");

    // a field constant in space and time is reproduced bitwise (nested linear interpolation)
    field C;
    for(double &c : C.c) c = 0.0;
    C.c[0] = 15.0;  C.c[4] = -2.5;  C.c[8] = 101325.0;
    write_file("met_test_tmp.r3m",7,5,0.0,0.0,100.0,100.0,times,C);
    met_file mc;
    mc.open("met_test_tmp.r3m",err);
    bool exact = true;
    for(double t : {0.0, 123.456, 777.7, 1199.9})
    {
        mc.advance(t,err);
        for(double x=-50.0; x<=650.0; x+=37.3)
        for(double y=-20.0; y<=420.0; y+=41.1)
        {
            mc.at(t,x,y,u,v,pa);
            exact = exact && u==15.0 && v==-2.5 && pa==101325.0;
        }
    }
    check(exact, "a field constant in space and time: bitwise the constant everywhere (as the uniform source)");

    // sources
    met_source_uniform s(20.0,90.0,100500.0);
    s.at(0.0,0.0,0.0,u,v,pa);
    check(std::fabs(u)<1.0e-14 && v==20.0 && pa==100500.0, "uniform source: 20 m/s to 90 deg = +y");
    met_source_uniform s0(20.0,0.0,101325.0);
    s0.at(0.0,0.0,0.0,u,v,pa);
    check(u==20.0 && v==0.0 && std::sqrt(u*u+v*v)==20.0, "uniform source along x: |U10| bitwise U10 (as A 571)");

    met_source_file sf;
    check(sf.open("met_test_tmp.r3m",err) && sf.first_time()==0.0, "file source");
    std::remove("met_test_tmp.r3m");
}

static void test_bad_files()
{
    std::cout<<"bad files"<<std::endl;
    std::string err;
    met_file m;
    bool bad = !m.open("no_such_file.r3m",err);
    check(bad, "missing file: "+err);

    {
    std::ofstream o("met_test_tmp.r3m",std::ios::binary);
    o<<"NOTAMETFILE.............................................";
    }
    met_file m2;
    bad = !m2.open("met_test_tmp.r3m",err);
    check(bad, "bad magic: "+err);

    field F;
    write_file("met_test_tmp.r3m",3,3,0.0,0.0,1.0,1.0,{0.0,10.0},F);
    {
    std::ofstream o("met_test_tmp.r3m",std::ios::binary|std::ios::app);
    o<<"x";
    }
    met_file m3;
    bad = !m3.open("met_test_tmp.r3m",err);
    check(bad, "wrong size: "+err);

    write_file("met_test_tmp.r3m",3,3,0.0,0.0,1.0,1.0,{0.0,10.0,10.0},F);
    met_file m4;
    bad = !m4.open("met_test_tmp.r3m",err);
    check(bad, "times not increasing: "+err);
    std::remove("met_test_tmp.r3m");
}

int main()
{
    test_datetime();
    test_drag();
    test_frame();
    test_python_file();
    test_window();
    test_uniform_and_moving();
    test_bad_files();

    std::cout<<std::endl<<(nfail==0 ? "met_test: all passed" : "met_test: "+std::to_string(nfail)+" FAILED")<<std::endl;
    return nfail==0 ? 0 : 1;
}
