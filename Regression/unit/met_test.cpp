// Architect: Hans Bihs
// Standalone verification of the REEF3D atmospheric forcing kernels (W 200): date-times, drag laws (met_drag.h),
// geographic frame (met_frame) against tools/seastate_forcing.py, binary met file reader (met_file): header,
// file written by tools/met_forcing.py (met_test.r3m), bilinear/linear interpolation exact for linear fields,
// windows (only the needed part read, the same values), record advance, clamping in space and time, uniform
// series, moving grids, bad files; uniform and file sources (met_source.h); the high-wind drag laws W 207 5, 6; the
// parametric vortex (met_vortex: Holland 1980, GAHM, Holland 2010) against tools/met_track.py (met_vortex_ref.txt,
// met_vortex_nh.trk, met_vortex_sh.trk), its profiles, the vortex source on the model grid; the blend W 200 4
// (met_blend: background, moving nest, vortex) against tools/met_forcing.py (met_blend_ref.txt, met_blend_bg.r3m,
// met_blend_nest.r3m, written by met_blend_make.py), the same values with rank windows, the layers bitwise at weights
// 0 and 1, the tapers, the ring pressure read once per record, the environment W 226 (background = vortex + linear field). No MPI, no REEF3D binary.
// Build:  g++ -O2 -std=c++20 -I../../src met_test.cpp ../../src/met_file.cpp ../../src/met_frame.cpp ../../src/met_vortex.cpp ../../src/met_blend.cpp -o met_test
// Run:    ./met_test
#include"met_drag.h"
#include"met_file.h"
#include"met_frame.h"
#include"met_source.h"
#include"met_vortex.h"
#include"met_blend.h"
#include<sstream>
#include<map>
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


static void test_drag_highwind()
{
    std::cout<<"high-wind drag laws (W 207 5, 6)"<<std::endl;
    bool ok = true;
    for(double U=0.5; U<90.0; U+=0.5)
    {
        const double g = 0.001*(0.75+0.067*U);
        ok = ok && met_drag(5,U)==std::min(g,2.5e-3) && met_drag(5,U,2.0e-3)==std::min(g,2.0e-3);
    }
    check(ok, "W 207 5: Garratt capped at W 213 (2.5e-3 and 2.0e-3), 0.5 - 90 m/s");
    check(std::fabs(met_drag(5,20.0)-2.09e-3)<1.0e-15 && met_drag(5,30.0)==2.5e-3, "W 207 5: below the cap Garratt, above the cap");
    double cmax = 0.0, umax = 0.0, cmin = 1.0;
    bool mono = true;
    for(double U=20.0; U<90.0; U+=0.1)
    {
        const double c = met_drag(6,U);
        if(c>cmax) {cmax = c; umax = U;}
        if(U>35.0) cmin = std::min(cmin,c);
        if(U<31.0) mono = mono && met_drag(6,U+0.1)>=c;
    }
    check(std::fabs(cmax-2.03e-3)<2.0e-6 && std::fabs(umax-31.4)<0.2, "W 207 6: peak 2.03e-3 at 31.4 m/s (Zijlema et al.)");
    check(std::fabs(cmin-1.5e-3)<1.0e-15 && met_drag(6,80.0)==1.5e-3 && met_drag(6,50.0)==met_drag(4,50.0), "W 207 6: not below the floor W 214 beyond the peak, Zijlema above it");
    check(mono && met_drag(6,10.0)==met_drag(4,10.0), "W 207 6: Zijlema et al. up to the peak");
    check(met_drag_law(0,1)==1 && met_drag_law(0,2)==1 && met_drag_law(0,3)==5 && met_drag_law(4,3)==4, "W 207 0: law 1 for W 200 1/2, 5 for the vortex, a given law kept");
}

static void test_vortex()
{
    std::cout<<"parametric vortex (W 200 3) against tools/met_track.py"<<std::endl;
    std::map<std::string,met_vortex> vx;
    std::string err;
    for(const char *f : {"met_vortex_nh.trk","met_vortex_sh.trk"})
    {
        met_vortex v;
        const bool ok = v.open(f,err);
        check(ok && v.rec.size()==5, std::string("open ")+f+(ok ? "" : ": "+err));
        vx[f] = v;
    }

    std::ifstream in("met_vortex_ref.txt");
    std::string line;
    int n = 0;
    double eu = 0.0, ep = 0.0;
    std::map<int,double> em;
    while(std::getline(in,line))
    {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream ls(line);
        std::string fn;
        int model, inflow, trans;
        double t, lon, lat, u, v, p;
        ls>>fn>>model>>inflow>>trans>>t>>lon>>lat>>u>>v>>p;
        met_vortex &V = vx[fn];
        V.par.model = model;  V.par.inflow = inflow;  V.par.trans = trans;
        V.par.wrf = 0.9;  V.par.one2ten = 0.8928;
        double uc, vc, pc;
        V.at(t,lon,lat,uc,vc,pc);
        const double e = std::max(std::fabs(uc-u),std::fabs(vc-v));
        eu = std::max(eu,e);
        ep = std::max(ep,std::fabs(pc-p));
        em[model] = std::max(em[model],e);
        ++n;
    }
    check(n==1200, "1200 reference points (2 tracks, NH and SH, models 1-3, inflow and translation options, 4 times)");
    char buf[160];
    snprintf(buf,sizeof(buf),"wind = Python to 1e-9 m/s (max %.2e m/s), pressure to 1e-7 Pa (max %.2e Pa)",eu,ep);
    check(eu<1.0e-9 && ep<1.0e-7, buf);

    // profiles
    const double f = 2.0*7.2921159e-5*std::sin(26.0*3.14159265358979323846/180.0);
    double vmx = 0.0, rmx = 0.0;
    for(double r=5000.0; r<200000.0; r+=10.0)
    {
        const double v = met_vortex::gahm(r,50.0,40000.0,1.6,1.03,f);
        if(v>vmx) {vmx = v; rmx = r;}
    }
    check(std::fabs(met_vortex::gahm(40000.0,50.0,40000.0,1.6,1.03,f)-50.0)<1.0e-12, "GAHM: Vg(Rmax) = Vmax exactly");
    check(std::fabs(rmx-40000.0)<2000.0, "GAHM: the maximum of the profile near Rmax ("+std::to_string(rmx)+" m)");
    check(std::fabs(met_vortex::holland10(30000.0,45.0,30000.0,1.5,0.8,200000.0)-45.0)<1.0e-12, "Holland 2010: V(Rmax) = Vms");
    check(met_vortex::fang(10.0,20.0)==5.0 && std::fabs(met_vortex::fang(22.0,20.0)-17.5)<1.0e-12 && met_vortex::fang(30.0,20.0)==25.0, "inflow angle of ADCIRC: 10 r/Rm, 10 + 75 (r/Rm - 1), 25 deg");
    const double B = 1.4, dp = 5000.0, rho = 1.15;
    const double vc = std::sqrt(B*dp/rho/std::exp(1.0));
    check(std::fabs(met_vortex::holland80(30000.0,30000.0,B,dp,rho,0.0)-vc)<1.0e-12, "Holland 1980 without f: Vmax = sqrt(B dp/(rho e)) at Rmax");

    // hemispheres: counter-clockwise in the north, clockwise in the south (no translation, no inflow)
    met_vortex &N = vx["met_vortex_nh.trk"], &S = vx["met_vortex_sh.trk"];
    N.par.inflow = 0; N.par.trans = 0; N.par.model = 2;
    S.par.inflow = 0; S.par.trans = 0; S.par.model = 2;
    double lon, lat, ue, vn, pa;
    const met_vortex_rec &rn = N.rec[0], &rs = S.rec[0];
    // point east of the centre: NH wind to the north, SH to the south
    N.eval(rn,rn.lon,rn.lat,0.0,0.0,rn.lon+0.5,rn.lat,ue,vn,pa);
    double ue2, vn2, pa2;
    S.eval(rs,rs.lon,rs.lat,0.0,0.0,rs.lon+0.5,rs.lat,ue2,vn2,pa2);
    check(vn>0.0 && std::fabs(ue)<1.0e-2*vn && vn2<0.0 && std::fabs(ue2)<1.0e-2*std::fabs(vn2), "east of the centre: wind to the north (NH), to the south (SH)");
    check(pa<rn.pn && pa>rn.pc, "pressure between pc and pn");
    N.eval(rn,rn.lon,rn.lat,0.0,0.0,rn.lon+0.005,rn.lat,ue,vn,pa);
    check(ue==0.0 && vn==0.0 && pa==rn.pc, "inside 1 km: no wind, pc");

    // the vortex source on a model grid (frame about the track): bitwise the geographic evaluation turned
    met_frame fr(-86.0,26.0,15.0,0.0,0.0);
    met_source_vortex src(fr);
    src.open("met_vortex_nh.trk",err);
    src.v.par = N.par;
    double x, y, u, v, p;
    fr.to_model(-85.3,26.4,x,y);
    src.at(rn.t+3600.0,x,y,u,v,p);
    fr.to_geo(x,y,lon,lat);
    src.v.at(rn.t+3600.0,lon,lat,ue,vn,pa);
    const double chi = fr.north_geo(lon,lat);
    check(std::fabs(std::hypot(u,v)-std::hypot(ue,vn))<1.0e-12 && p==pa, "vortex source: speed and pressure kept by the turn to the model axes");
    check(std::fabs(std::atan2(v,u) - (std::atan2(vn,ue) + chi - 0.5*3.14159265358979323846))<1.0e-9
          || std::fabs(std::fabs(std::atan2(v,u) - (std::atan2(vn,ue) + chi - 0.5*3.14159265358979323846)) - 2.0*3.14159265358979323846)<1.0e-9,
          "vortex source: direction turned by the angle of true north in the model frame (rotation 15 deg)");

    // bad tracks
    {
    std::ofstream o("met_test_tmp.trk");
    o<<"REEF3D-MET-TRACK 1 2 1.15\n20181010.000000 1 2 3\n";
    }
    met_vortex bad;
    const bool b1 = !bad.open("met_test_tmp.trk",err);
    check(b1, "incomplete record refused: "+err);
    std::remove("met_test_tmp.trk");
}

// ---------------------------------------------------------------------------------------------
// blend W 200 4 (met_blend) against tools/met_forcing.py Composite
// ---------------------------------------------------------------------------------------------
static met_blend_par blend_variant(char v, met_vortex_par &vp)
{
    met_blend_par b;
    vp = met_vortex_par();
    b.nest = 1;
    b.vortex = 1;

    if(v=='B')
    {
    b.shape = 2;  b.pressure = 0;  b.basis = 2;  b.r1 = 0.5;  b.r2 = 1.0;  b.vortex_time = 0.0;
    b.nest_width = 30.0e3;  b.nest_time = 0.0;  b.nest_relative = 1;
    vp.model = 1;  vp.trans = 2;
    }

    if(v=='C')
    {
    b.nest = 0;  b.basis = 3;  b.r1 = 80.0e3;  b.r2 = 250.0e3;  b.env = 1;
    vp.model = 3;  vp.inflow = 0;
    }

    if(v=='D')
    {
    b.vortex = 0;  b.nest_relative = 1;  b.nest_width = 20.0e3;
    }

    return b;
}

static void test_blend()
{
    std::cout<<"blend W 200 4 (met_blend) against tools/met_forcing.py (met_blend_ref.txt)"<<std::endl;

    const met_frame fr(-86.0,26.0,10.0,5000.0,-3000.0);
    std::string err;

    struct row {char v; double t, x, y, u, v_, pa, wn, wv;};
    std::vector<row> ref;
    {
    std::ifstream in("met_blend_ref.txt");
    std::string line;

        while(std::getline(in,line))
        {
            if(line.empty() || line[0]=='#')
            continue;

        std::istringstream ls(line);
        row r;
        ls>>r.v>>r.t>>r.x>>r.y>>r.u>>r.v_>>r.pa>>r.wn>>r.wv;
        ref.push_back(r);
        }
    }
    check(ref.size()==1920, "reference: "+std::to_string(ref.size())+" rows (120 points, 4 times, 4 parameter sets)");

    double du = 0.0, dpa = 0.0;
    int n = 0, nbit = 0;

    for(char var : {'A','B','C','D'})
    {
    met_vortex_par vp;
    met_source_blend src(fr,blend_variant(var,vp));
    const bool ok = src.open("met_blend_bg.r3m","met_blend_nest.r3m","met_vortex_nh.trk",err);
    src.vx.par = vp;

        if(!ok)
        {
        check(false, "open: "+err);
        return;
        }

    double tlast = -1.0;

        for(const row &r : ref)
        {
            if(r.v!=var)
            continue;

            if(r.t!=tlast)
            {
            src.advance(r.t,err);
            tlast = r.t;
            }

        double u, v, pa;
        src.at(r.t,r.x,r.y,u,v,pa);

        du = std::max(du, std::max(std::fabs(u-r.u),std::fabs(v-r.v_)));
        dpa = std::max(dpa, std::fabs(pa-r.pa));
        nbit += (u==r.u && v==r.v_ && pa==r.pa);
        ++n;
        }
    }
    char buf[256];
    snprintf(buf,sizeof(buf),"wind = Python to 1e-9 m/s (max %.2e m/s), pressure to 1e-6 Pa (max %.2e Pa); %d of %d bitwise",du,dpa,nbit,n);
    check(n==1920 && du<1.0e-9 && dpa<1.0e-6, buf);

    // rank windows: a source that reads only a box gives the same values in the box (ring from the file,
    // the nest window grown by the move of the storm between two nest records for W 227 1)
    for(char var : {'A','D'})
    {
    met_vortex_par vp;
    met_source_blend all(fr,blend_variant(var,vp)), part(fr,blend_variant(var,vp));
    all.vx.par = vp;
    part.vx.par = vp;
    all.open("met_blend_bg.r3m","met_blend_nest.r3m","met_vortex_nh.trk",err);
    part.open("met_blend_bg.r3m","met_blend_nest.r3m","met_vortex_nh.trk",err);

    bool same = true;
    const double t = met_datetime(20181010.0915);
    all.advance(t,err);

    // A: a large box; D (nest only): a box of 100 x 50 km on the storm (it moves 60 km between the nest records)
    double lonc, latc, xs, ys;
    all.vx.centre(t,lonc,latc);
    fr.to_model(lonc,latc,xs,ys);
    const double bx0 = var=='A' ? -100.0e3 : xs - 50.0e3, bx1 = var=='A' ? 150.0e3 : xs + 50.0e3;
    const double by0 = var=='A' ? 100.0e3 : ys - 25.0e3, by1 = var=='A' ? 400.0e3 : ys + 25.0e3;
    part.window(bx0,bx1,by0,by1);
    part.advance(t,err);

        for(double x=bx0; x<=bx1; x+=(bx1-bx0)/20.0)
        for(double y=by0; y<=by1; y+=(by1-by0)/20.0)
        {
        double u1, v1, p1, u2, v2, p2;
        all.at(t,x,y,u1,v1,p1);
        part.at(t,x,y,u2,v2,p2);
        same = same && u1==u2 && v1==v2 && p1==p2;
        }
    check(same && all.env.p==part.env.p && all.env.u==part.env.u && all.env.gx==part.env.gx,
          std::string("parameter set ")+var+(var=='A' ? ": a rank window gives bitwise the same composite (ring read from the file, not the window)"
                                                      : ": the same with the storm-relative nest (W 227 1, nest window grown by the move of the storm)"));
    check(part.bg.values_read < all.bg.values_read/4, "the window reads less of the background ("+std::to_string(part.bg.values_read)+" of "
                                                     +std::to_string(all.bg.values_read)+" values)");

        if(var=='D')
        continue;

    // the ring pressure is read once per background record
    const long long s0 = all.ring_samples;
        for(int k=1; k<=20; ++k)
        all.advance(t + 60.0*k,err);
    check(all.ring_samples==s0, "ring pressure: no new reads within the same background records");
    all.advance(met_datetime(20181010.1430),err);
    check(all.ring_samples==s0+1, "ring pressure: one new record read when the background moves on");
    }

    // layers bitwise: W 220 0 0 = the background file; far beyond r2 = the background; inside r1 = the vortex
    {
    met_blend_par b;
    b.nest = 0;  b.vortex = 0;
    met_source_blend bgo(fr,b);
    met_source_file fs;
    bgo.open("met_blend_bg.r3m","","",err);
    fs.open("met_blend_bg.r3m",err);
    bool same = true;

        for(double t : {met_datetime(20181009.2), met_datetime(20181010.07), met_datetime(20181010.2)})
        {
        bgo.advance(t,err);
        fs.advance(t,err);

            for(double x=-800.0e3; x<=800.0e3; x+=37.0e3)
            for(double y=-600.0e3; y<=1000.0e3; y+=41.0e3)
            {
            double u1, v1, p1, u2, v2, p2;
            bgo.at(t,x,y,u1,v1,p1);
            fs.at(t,x,y,u2,v2,p2);
            same = same && u1==u2 && v1==v2 && p1==p2;
            }
        }
    check(same, "W 220 0 0: bitwise the met file source (W 200 2)");

    met_vortex_par vp;
    met_blend_par bv;
    bv.vortex = 1;  bv.basis = 3;  bv.r1 = 1.0e8;  bv.r2 = 2.0e8;  bv.pressure = 0;  bv.vortex_time = 0.0;
    met_source_blend vo(fr,bv);
    met_source_vortex sv(fr);
    vo.open("met_blend_bg.r3m","","met_vortex_nh.trk",err);
    sv.open("met_vortex_nh.trk",err);
    same = true;

        for(double t : {met_datetime(20181010.0), met_datetime(20181010.0730), met_datetime(20181011.0)})
        {
        vo.advance(t,err);

            for(double x=-500.0e3; x<=700.0e3; x+=53.0e3)
            for(double y=-300.0e3; y<=900.0e3; y+=47.0e3)
            {
            double u1, v1, p1, u2, v2, p2;
            vo.at(t,x,y,u1,v1,p1);
            sv.at(t,x,y,u2,v2,p2);
            same = same && u1==u2 && v1==v2 && p1==p2;
            }
        }
    check(same, "weight 1 everywhere (r1 = 1e8 m, W 223 0, W 224 0): bitwise the vortex source (W 200 3)");

    vo.advance(met_datetime(20181009.2),err);
    double u1, v1, p1, u2, v2, p2;
    vo.at(met_datetime(20181009.2),0.0,0.0,u1,v1,p1);
    met_source_file fb;
    fb.open("met_blend_bg.r3m",err);
    fb.advance(met_datetime(20181009.2),err);
    fb.at(met_datetime(20181009.2),0.0,0.0,u2,v2,p2);
    check(u1==u2 && v1==v2 && p1==p2 && vo.wv_t==0.0, "before the first track record: no vortex, the background");
    }

    // tapers
    {
    const double e = 1.0e-6;
    const double d0 = (met_source_blend::taper(e,1) - met_source_blend::taper(0.0,1))/e;
    const double d1 = (met_source_blend::taper(1.0,1) - met_source_blend::taper(1.0-e,1))/e;
    bool mono = true;
        for(int k=1; k<=100; ++k)
        mono = mono && met_source_blend::taper(0.01*k,1)>=met_source_blend::taper(0.01*(k-1),1);
    check(met_source_blend::taper(-0.3,1)==0.0 && met_source_blend::taper(1.7,1)==1.0 && std::fabs(met_source_blend::taper(0.5,1)-0.5)<1.0e-15
          && d0<1.0e-5 && d1<1.0e-5 && mono, "cosine taper: 0 to 1, monotonic, zero slope at both ends (W 225 1)");
    check(met_source_blend::taper(0.3,2)==0.3 && met_source_blend::taper(1.2,2)==1.0, "linear taper (W 225 2, ADCIRC NWS 29/30)");
    check(met_source_blend::taper_time(5.0,10.0,20.0,3.0,1)==0.0 && met_source_blend::taper_time(15.0,10.0,20.0,3.0,1)==1.0
          && met_source_blend::taper_time(11.5,10.0,20.0,3.0,1)==0.5*(1.0-std::cos(pi*0.5)) && met_source_blend::taper_time(20.0,10.0,20.0,0.0,1)==1.0
          && met_source_blend::taper_time(20.5,10.0,20.0,0.0,1)==0.0, "taper in time: 0 outside the records, over T after the first and before the last");
    }

    // the vortex pressure rescaled to the ring: pv' = pb_ring at the ring where the deficit of the vortex is ~0
    {
    met_vortex_par vp;
    met_source_blend a(fr,blend_variant('A',vp));
    a.open("met_blend_bg.r3m","met_blend_nest.r3m","met_vortex_nh.trk",err);
    const double t = met_datetime(20181010.12);
    a.advance(t,err);
    double sum = 0.0;
        for(int m=0; m<16; ++m)
        {
        double lon2, lat2, x, y, val;
        met_vortex::destination(a.lonc,a.latc,22.5*m,a.r2m,lon2,lat2);
        fr.to_model(lon2,lat2,x,y);
        a.bg.sample(a.bg.nrec>1 ? 3 : 0,2,x,y,val,err);
        sum += val;
        }
    check(std::fabs(a.env.p - sum/16.0)<1.0e-9 && a.env.p<a.pn, "ring pressure: the mean background pressure on the ring r2 at the record time ("
          +std::to_string(a.env.p)+" Pa, pn of the track "+std::to_string(a.pn)+" Pa)");
    check(std::fabs(a.r2m - 6.0*a.rmw)<1.0e-6 && std::fabs(a.r1m - 3.0*a.rmw)<1.0e-6, "radii W 222 1 3 6: 3 and 6 rmw");
    }

    // environment W 226 1: background = the vortex + a field linear in x, y, t -> the environment is the linear field
    {
    field F;
    const double t0 = met_datetime(20181010.0);
    const std::vector<double> times = {t0, t0+21600.0};
    met_source_vortex sv(fr);
    sv.open("met_vortex_nh.trk",err);
    {
    std::ofstream out("met_test_env.r3m", std::ios::binary);
    out.write("REEF3DMT",8);
    const int nx = 121, ny = 121;
    const double x0 = -700.0e3, y0 = -500.0e3, dx = 10.0e3;
    int32_t iv[8] = {1, nx, ny, 3, 0, 2, 0, 0};
    out.write((char*)iv,sizeof(iv));
    double dv[13] = {x0, y0, dx, dx, 0.0, 0.0, 0.0, 0.0, 0.0, 600.0, 0.0, 0.0, 0.0};
    out.write((char*)dv,sizeof(dv));
        for(double tn : times)
        {
        double h[4] = {tn, x0, y0, 0.0};
        out.write((char*)h,sizeof(h));
        std::vector<float> a[3];
            for(auto &q : a)
            q.resize(size_t(nx)*ny);
            for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
            {
            double u, v, pa, uv, vv, pv;
            F.at(x0 + i*dx, y0 + j*dx, tn - t0, u, v, pa);
            sv.at(tn, x0 + i*dx, y0 + j*dx, uv, vv, pv);
            a[0][j*nx+i] = float(u + uv);
            a[1][j*nx+i] = float(v + vv);
            a[2][j*nx+i] = float(pa + pv - 101200.0);
            }
            for(auto &q : a)
            out.write((char*)q.data(), 4*q.size());
        }
    }
    met_blend_par b;
    b.env = 1;
    b.vortex_time = 0.0;
    met_source_blend e(fr,b);
    e.open("met_test_env.r3m","","met_vortex_nh.trk",err);
    const double t = t0 + 7200.0;
    e.advance(t,err);

    // the linear field at the mean of the ring points, linear in time as the ring
    double mxy[2][2];
        for(int q=0; q<2; ++q)
        {
        double lon, lat, p0, p1, rm, r3;
        e.vx.state(times[q],lon,lat,p0,p1,rm,r3);
        mxy[q][0] = mxy[q][1] = 0.0;
            for(int m=0; m<16; ++m)
            {
            double lon2, lat2, x, y;
            met_vortex::destination(lon,lat,22.5*m,6.0*rm,lon2,lat2);
            fr.to_model(lon2,lat2,x,y);
            mxy[q][0] += x/16.0;
            mxy[q][1] += y/16.0;
            }
        }
    double u0, v0, p0, u1, v1, p1;
    F.at(mxy[0][0],mxy[0][1],0.0,u0,v0,p0);
    F.at(mxy[1][0],mxy[1][1],21600.0,u1,v1,p1);
    const double w = 7200.0/21600.0;
    const double ue = u0+w*(u1-u0), ve = v0+w*(v1-v0);
    char buf[300];
    snprintf(buf,sizeof(buf),"W 226 1, background = vortex + linear field: environment gradient %.6f %.6f Pa/m (linear -0.01 0.02), wind %.3f %.3f m/s (linear %.3f %.3f)",
             e.env.gx,e.env.gy,e.env.u,e.env.v,ue,ve);
    check(std::fabs(e.env.gx+0.01)<2.0e-4 && std::fabs(e.env.gy-0.02)<2.0e-4 && std::fabs(e.env.u-ue)<0.05 && std::fabs(e.env.v-ve)<0.05, buf);

    // inside r1 the composite is the vortex plus the environment
    double x, y, uc, vc, pcmp;
    fr.to_model(e.lonc+0.4,e.latc+0.2,x,y);
    e.at(t,x,y,uc,vc,pcmp);
    met_blend_par b0 = b;
    b0.env = 0;
    met_source_blend e0(fr,b0);
    e0.open("met_test_env.r3m","","met_vortex_nh.trk",err);
    e0.advance(t,err);
    double ud, vd, pd;
    e0.at(t,x,y,ud,vd,pd);
    check(std::fabs(uc-ud-e.env.u)<1.0e-9 && std::fabs(vc-vd-e.env.v)<1.0e-9
          && std::fabs(pcmp-pd-(e.env.gx*(x-e.xcm)+e.env.gy*(y-e.ycm)))<1.0e-6, "W 226 1 inside r1: the vortex plus the environment wind and pressure plane");
    check(std::fabs(e.env.p - e.env.pv - (p0+w*(p1-p0) - 101200.0))<5.0, "W 223 1: background minus vortex on the ring = the linear field ("
          +std::to_string(e.env.p - e.env.pv)+" Pa)");
    std::remove("met_test_env.r3m");
    }

    // bad input
    {
    met_blend_par b;
    b.nest = 1;
    met_source_blend x(fr,b);
    const bool refused = !x.open("met_blend_bg.r3m","no_such_nest.r3m","met_vortex_nh.trk",err);
    check(refused, "missing nest file refused: "+err);
    }
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
    test_drag_highwind();
    test_vortex();
    test_blend();

    std::cout<<std::endl<<(nfail==0 ? "met_test: all passed" : "met_test: "+std::to_string(nfail)+" FAILED")<<std::endl;
    return nfail==0 ? 0 : 1;
}
