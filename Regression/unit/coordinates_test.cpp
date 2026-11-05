// Standalone check of the world <-> model coordinate transform (class coordinates, used when the
// grid is rotated/shifted, cms_flag==1): no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src coordinates_test.cpp -o coordinates_test
// Run:    ./coordinates_test
// Checks: XYout(XYin(x,y)) == (x,y) and XYin(XYout(x,y)) == (x,y) over a range of grid angles and
// origins, XYin agrees with Xin/Yin and XYout with Xout/Yout, a hand-computed rotation by 90 deg,
// lengths are kept, angles in/out are inverse, and cms_flag==0 leaves everything unchanged.
// Architect: Hans Bihs

// coordinates.cpp needs only these members of the lexer and nothing of increment: small
// stand-ins, so the test builds without the solver (the include guards keep the real headers out)
#include<cmath>
#define LEXER_H_
#define INCREMENT_H_
inline constexpr double PI = 3.14159265359;
class increment
{
public:
    virtual ~increment() {}
};
class lexer
{
public:
    int cms_flag = 0;
    double global_orig_x = 0.0, global_orig_y = 0.0, alpha_grid = 0.0;
};
#include"../../src/coordinates.cpp"

#include<cstdio>
#include<iostream>
#include<string>

static int nfail = 0, npass = 0;
static void check(bool ok, const std::string &what)
{
    if(ok) ++npass;
    else
    {
        ++nfail;
        if(nfail<=20)
        std::cout<<"  FAIL  "<<what<<std::endl;
        if(nfail==20)
        std::cout<<"  (further failures not printed)"<<std::endl;
    }
}

static bool near(double a, double b, double tol=1.0e-9)
{
    return std::fabs(a-b) <= tol*(1.0+std::fabs(a)+std::fabs(b));
}

static std::string tag(double alpha, double x0, double y0, double x, double y)
{
    char buf[160];
    std::snprintf(buf,sizeof(buf),"alpha=%g deg, orig=(%g,%g), point=(%g,%g)",alpha*180.0/PI,x0,y0,x,y);
    return buf;
}

int main()
{
    lexer lex;
    coordinates c(&lex);

    // cms_flag==0: identity
    {
        lex.cms_flag = 0;
        lex.global_orig_x = 5.0; lex.global_orig_y = -3.0; lex.alpha_grid = 0.7;
        double x = 12.5, y = -4.25;
        c.XYin(x,y);
        check(x==12.5 && y==-4.25, "cms_flag=0: XYin leaves the point unchanged");
        c.XYout(x,y);
        check(x==12.5 && y==-4.25, "cms_flag=0: XYout leaves the point unchanged");
        check(c.Xin(1.0,2.0)==1.0 && c.Yin(1.0,2.0)==2.0, "cms_flag=0: Xin/Yin unchanged");
        check(c.Alpha_deg_in(30.0)==30.0, "cms_flag=0: Alpha_deg_in unchanged");
    }

    lex.cms_flag = 1;

    // hand-computed: grid rotated by +90 deg, origin (10,20): world (10,21) is model (1,0)
    {
        lex.global_orig_x = 10.0; lex.global_orig_y = 20.0; lex.alpha_grid = 0.5*PI;
        double x = 10.0, y = 21.0;
        c.XYin(x,y);
        check(near(x,1.0) && std::fabs(y)<1.0e-9, "90 deg: XYin(10,21) == (1,0)");
        x = 1.0; y = 0.0;
        c.XYout(x,y);
        check(near(x,10.0) && near(y,21.0), "90 deg: XYout(1,0) == (10,21)");
    }

    // pure shift: alpha=0
    {
        lex.global_orig_x = 100.0; lex.global_orig_y = -50.0; lex.alpha_grid = 0.0;
        double x = 103.0, y = -46.0;
        c.XYin(x,y);
        check(near(x,3.0) && near(y,4.0), "shift only: XYin subtracts the origin");
    }

    // round trips and consistency over a sweep of angles, origins and points
    const double alphas[] = {0.0, 1.0e-3, 0.3, 0.25*PI, 0.5*PI, 2.0, PI, -0.6, -0.5*PI, 5.5};
    const double origins[][2] = {{0.0,0.0}, {10.0,20.0}, {-347.25,1250.5}, {4.6e5,7.03e6}};
    const double points[][2] = {{0.0,0.0}, {1.0,0.0}, {0.0,1.0}, {12.5,-7.75}, {-300.0,42.0}, {4.6e5+812.0,7.03e6-95.5}};

    int ncase = 0;
    for(double a : alphas)
    for(auto &o : origins)
    for(auto &pt : points)
    {
        lex.alpha_grid = a; lex.global_orig_x = o[0]; lex.global_orig_y = o[1];
        const double xw = pt[0], yw = pt[1];
        const std::string t = tag(a,o[0],o[1],xw,yw);
        const double scale = 1.0e-12*(std::fabs(o[0])+std::fabs(o[1])+std::fabs(xw)+std::fabs(yw)+1.0);

        // world -> model -> world
        double x = xw, y = yw;
        c.XYin(x,y);
        const double xm = x, ym = y;
        c.XYout(x,y);
        check(std::fabs(x-xw)<=scale && std::fabs(y-yw)<=scale, "XYout(XYin(p)) == p, "+t);

        // model -> world -> model (treat the world point as a model point)
        x = xw; y = yw;
        c.XYout(x,y);
        c.XYin(x,y);
        check(std::fabs(x-xw)<=scale && std::fabs(y-yw)<=scale, "XYin(XYout(p)) == p, "+t);

        // XYin agrees with Xin/Yin, XYout with Xout/Yout
        check(std::fabs(xm-c.Xin(xw,yw))<=scale && std::fabs(ym-c.Yin(xw,yw))<=scale, "XYin == (Xin,Yin), "+t);
        x = xw; y = yw;
        c.XYout(x,y);
        check(std::fabs(x-c.Xout(xw,yw))<=scale && std::fabs(y-c.Yout(xw,yw))<=scale, "XYout == (Xout,Yout), "+t);

        // a rigid transform keeps the distance to the origin
        const double rw = std::hypot(xw-o[0], yw-o[1]);
        check(std::fabs(std::hypot(xm,ym)-rw)<=scale, "|XYin(p)| == |p-orig|, "+t);
        ++ncase;
    }

    // a segment keeps its direction relative to the grid: direction in = direction out - alpha
    {
        lex.alpha_grid = 0.4; lex.global_orig_x = 7.0; lex.global_orig_y = -2.0;
        double xs = 1.0, ys = 1.0, xe = 4.0, ye = 5.0;
        const double dir_world = std::atan2(ye-ys, xe-xs);
        c.XYin(xs,ys); c.XYin(xe,ye);
        const double dir_model = std::atan2(ye-ys, xe-xs);
        check(near(dir_model, c.Alpha_rad_in(dir_world)), "segment direction: model == Alpha_rad_in(world)");
        check(near(c.Alpha_rad_out(c.Alpha_rad_in(dir_world)), dir_world), "Alpha_rad_out(Alpha_rad_in(a)) == a");
        check(near(c.Alpha_deg_out(c.Alpha_deg_in(33.0)), 33.0), "Alpha_deg_out(Alpha_deg_in(a)) == a");
    }

    std::cout<<"coordinates_test: "<<ncase<<" transform cases, "<<npass<<" checks passed, "<<nfail<<" failed"<<std::endl;
    return nfail==0 ? 0 : 1;
}
