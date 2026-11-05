// Standalone check of the prescribed body deformation of a 6DOF body (X 195, src/6DOF_deform.cpp):
// no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY deform_test.cpp -o deform_test
// Run:    ./deform_test
// Checks: mass distribution along x from the closed hull (divergence theorem) for a box and for the
// 2D fish profile of Kern & Koumoutsakos (2006) extruded in y (volume, bin volumes against the
// width profile, also with the hull oriented inwards); head, tail and length; the envelopes
// (tail amplitudes of types 1-3); the deformation velocity is the time derivative of the
// displacement (with and without the ramp); recoil: no net lateral momentum and no angular momentum
// of the deformation (displacement and velocity) about the CoG, and without recoil the raw wave.
// Architect: Hans Bihs

// 6DOF_deform.cpp needs only the body-frame triangles of sixdof_geometry
#include<cmath>
#define INCREMENT_H_
class increment
{
public:
    virtual ~increment() {}
};
#include"../../src/6DOF_geometry.h"
sixdof_geometry::sixdof_geometry(int n) : id(n), tri_x(nullptr), tri_y(nullptr), tri_z(nullptr),
                                          tri_x0(nullptr), tri_y0(nullptr), tri_z0(nullptr), tricount(0) {}
sixdof_geometry::~sixdof_geometry() {}
#include"../../src/6DOF_deform.cpp"

#include<cstdio>
#include<iostream>
#include<string>
#include<vector>
#include<array>

static int nfail = 0, npass = 0;
static void check(bool ok, const std::string &what)
{
    if(ok) ++npass;
    else
    {
        ++nfail;
        if(nfail<=20)
        std::cout<<"  FAIL  "<<what<<std::endl;
    }
}

typedef std::array<double,3> P3;

struct hull
{
    std::vector<std::array<P3,3> > t;

    void tri(const P3 &a, const P3 &b, const P3 &c) {t.push_back({a,b,c});}
    void quad(const P3 &a, const P3 &b, const P3 &c, const P3 &d) {tri(a,b,c); tri(a,c,d);}

    // body-frame triangles of a geometry (shifted by the centroid in x as sixdof_geometry stores them)
    void to(sixdof_geometry &g, double shift, bool flip)
    {
        g.tricount = int(t.size());
        g.tri_x0 = new double*[g.tricount];
        g.tri_y0 = new double*[g.tricount];
        g.tri_z0 = new double*[g.tricount];
        for(int n=0; n<g.tricount; ++n)
        {
            g.tri_x0[n] = new double[3];
            g.tri_y0[n] = new double[3];
            g.tri_z0[n] = new double[3];
            for(int q=0; q<3; ++q)
            {
                const int qq = flip ? 2-q : q;
                g.tri_x0[n][q] = t[n][qq][0] - shift;
                g.tri_y0[n][q] = t[n][qq][1];
                g.tri_z0[n][q] = t[n][qq][2];
            }
        }
    }
};

// box [x0,x1] x [y0,y1] x [z0,z1], outward
static hull box(double x0, double x1, double y0, double y1, double z0, double z1)
{
    hull h;
    const P3 a{x0,y0,z0}, b{x1,y0,z0}, c{x1,y1,z0}, d{x0,y1,z0};
    const P3 e{x0,y0,z1}, f{x1,y0,z1}, g{x1,y1,z1}, k{x0,y1,z1};
    h.quad(a,d,c,b);    // z0, normal -z
    h.quad(e,f,g,k);    // z1, +z
    h.quad(a,b,f,e);    // y0, -y
    h.quad(d,k,g,c);    // y1, +y
    h.quad(a,e,k,d);    // x0, -x
    h.quad(b,c,g,f);    // x1, +x
    return h;
}

// half width of the 2D fish of Kern & Koumoutsakos (2006) (Gazzola et al. 2011, linear decrease
// in the middle part), s from the head, L = 1
static double fish_w(double s)
{
    const double wh = 0.04, sb = 0.04, st = 0.95, wt = 0.01;
    if(s<0.0 || s>1.0) return 0.0;
    if(s<sb) return std::sqrt(std::max(0.0, 2.0*wh*s - s*s));
    if(s<st) return wh - (wh - wt)*(s - sb)/(st - sb);
    return wt*(1.0 - s)/(1.0 - st);
}

// the profile in the x-z plane (head at x = 1, tail at x = 0, lateral z), extruded over y in [0,b]
static hull fish2D(int N, double b, std::vector<double> &xs)
{
    hull h;
    xs.clear();
    // cosine clustering at both ends
    for(int i=0; i<=N; ++i)
    xs.push_back(0.5*(1.0 - std::cos(M_PI*i/double(N))));

    for(int i=0; i<N; ++i)
    {
        const double xa = xs[i], xb = xs[i+1];
        const double wa = fish_w(1.0 - xa), wb = fish_w(1.0 - xb);
        // upper surface z = +w (outward +z), lower z = -w
        h.quad(P3{xa,0.0,wa}, P3{xb,0.0,wb}, P3{xb,b,wb}, P3{xa,b,wa});
        h.quad(P3{xa,0.0,-wa}, P3{xa,b,-wa}, P3{xb,b,-wb}, P3{xb,0.0,-wb});
        // end caps y = 0 (outward -y) and y = b (+y)
        h.quad(P3{xa,0.0,-wa}, P3{xb,0.0,-wb}, P3{xb,0.0,wb}, P3{xa,0.0,wa});
        h.quad(P3{xa,b,-wa}, P3{xa,b,wa}, P3{xb,b,wb}, P3{xb,b,-wb});
    }
    return h;
}

int main()
{
    std::cout<<"deform_test: prescribed body deformation (X 195)"<<std::endl;

    // ---- box: volume, uniform mass distribution, head/tail/length
    {
        hull hb = box(-0.3,0.5, -0.1,0.1, -0.05,0.05);
        for(int flip=0; flip<2; ++flip)
        {
            sixdof_geometry g(0);
            hb.to(g,0.0,flip==1);
            sixdof_deform d(1, 0.1, 1.0, 1.0, 0.0, 1, true);
            d.setup(g);
            const double V = 0.8*0.2*0.1;
            check(std::fabs(d.volume()-V)<1.0e-12*V, "box volume "+std::to_string(d.volume())+(flip?" (inward hull)":""));
            check(std::fabs(d.length()-0.8)<1.0e-14 && std::fabs(d.xhead()-0.5)<1.0e-14 && std::fabs(d.xtail()+0.3)<1.0e-14, "box length, head and tail");
            // the recoil of a uniform body: the fitted line through the uniform mass distribution
            d.set_time(0.37);
            double m0 = 0.0, m1 = 0.0, href = 0.0;
            const int M = 4000;
            for(int k=0; k<M; ++k)
            {
                const double x = -0.3 + (k+0.5)*0.8/M;
                m0 += d.disp(x);
                m1 += x*d.disp(x);
                href += std::fabs(d.disp(x));
            }
            check(std::fabs(m0)<1.0e-4*href && std::fabs(m1)<1.0e-4*href, "box: recoil, sum h = "+std::to_string(m0/href)+", sum x h = "+std::to_string(m1/href));
        }
    }

    // ---- 2D fish (extruded): volume and mass distribution against the width profile
    {
        std::vector<double> xs;
        const double b = 0.25;
        hull hf = fish2D(400,b,xs);

        // area of the profile (trapezoidal on the same points, exact for the polygon)
        double A = 0.0, Ax = 0.0;
        for(size_t i=0; i+1<xs.size(); ++i)
        {
            const double wa = fish_w(1.0-xs[i]), wb = fish_w(1.0-xs[i+1]);
            const double dx = xs[i+1]-xs[i];
            A  += dx*(wa + wb);
            Ax += dx*(wa + wb)*0.5*(xs[i]+xs[i+1]);
        }
        const double xc = Ax/A;

        for(int flip=0; flip<2; ++flip)
        {
            sixdof_geometry g(0);
            hf.to(g,xc,flip==1);
            sixdof_deform d(1, 0.125, 1.0, 1.0, 1.0, 2, true);
            d.setup(g);
            check(std::fabs(d.volume()-A*b)<1.0e-10, "fish volume "+std::to_string(d.volume())+" vs "+std::to_string(A*b));
            check(std::fabs(d.length()-1.0)<1.0e-12, "fish length");

            // recoil at several times, displacement and velocity
            for(double t : {0.1, 0.37, 0.8, 1.3, 2.71})
            {
                d.set_time(t);
                // integrals with the exact width (fine midpoint rule along x)
                double I0 = 0.0, I1 = 0.0, J0 = 0.0, J1 = 0.0, href = 0.0, vref = 0.0;
                const int M = 20000;
                for(int k=0; k<M; ++k)
                {
                    const double xx = (k+0.5)/double(M);
                    const double m = 2.0*fish_w(1.0-xx)*b/M;
                    const double x = xx - xc;
                    I0 += m*d.disp(x);
                    I1 += m*x*d.disp(x);
                    J0 += m*d.vel(x);
                    J1 += m*x*d.vel(x);
                    href += m*std::fabs(d.disp(x));
                    vref += m*std::fabs(d.vel(x));
                }
                check(std::fabs(I0)<2.0e-3*href && std::fabs(I1)<2.0e-3*href, "fish recoil, displacement, t = "+std::to_string(t)+": "+std::to_string(I0/href)+" "+std::to_string(I1/href));
                check(std::fabs(J0)<2.0e-3*vref && std::fabs(J1)<2.0e-3*vref, "fish recoil, velocity, t = "+std::to_string(t)+": "+std::to_string(J0/vref)+" "+std::to_string(J1/vref));
            }

            // the velocity is the time derivative of the displacement (in the ramp and after it)
            for(double t : {0.3, 0.7, 1.6, 3.2})
            {
                double emax = 0.0, vmax = 0.0;
                const double dt = 1.0e-6;
                for(int k=0; k<=50; ++k)
                {
                    const double x = d.xtail() + k*d.length()/50.0;
                    d.set_time(t+dt);
                    const double hp = d.disp(x);
                    d.set_time(t-dt);
                    const double hm = d.disp(x);
                    d.set_time(t);
                    emax = std::max(emax, std::fabs((hp-hm)/(2.0*dt) - d.vel(x)));
                    vmax = std::max(vmax, std::fabs(d.vel(x)));
                }
                check(emax<1.0e-6*vmax, "velocity = dh/dt at t = "+std::to_string(t)+": "+std::to_string(emax/vmax));
            }
        }

        // without recoil: the raw travelling wave (head at s = 0, tail amplitude 0.125 L)
        sixdof_geometry g(0);
        hf.to(g,xc,false);
        sixdof_deform d(1, 0.125, 1.0, 1.0, 0.0, 2, false);
        d.setup(g);
        double hmax = 0.0;
        for(int k=0; k<1000; ++k)
        {
            d.set_time(k/1000.0);
            hmax = std::max(hmax, std::fabs(d.disp(d.xtail())));
            const double s = 0.3;
            const double href = 0.125*(s + 0.03125)/1.03125*std::sin(2.0*M_PI*(s - k/1000.0));
            check(std::fabs(d.disp(d.xhead() - s) - href)<1.0e-14, "raw wave, type 1");
        }
        check(std::fabs(hmax-0.125)<1.0e-5, "tail amplitude type 1: "+std::to_string(hmax));
        d.set_time(0.0);
        check(std::fabs(d.disp(d.xhead()))<1.0e-15, "head at s = 0, t = 0");
    }

    // ---- envelopes: tail amplitude (r = 1) of the three types
    {
        hull hb = box(0.0,2.0, -0.1,0.1, -0.1,0.1);
        sixdof_geometry g(0);
        hb.to(g,1.0,false);
        sixdof_deform d1(1, 0.1, 1.0, 1.0, 0.0, 1, false), d2(2, 0.1, 1.0, 1.0, 0.0, 1, false), d3(3, 0.1, 1.0, 1.0, 0.0, 1, false);
        d1.setup(g); d2.setup(g); d3.setup(g);
        check(std::fabs(d1.tail_amplitude()-0.2)<1.0e-14, "type 1 tail amplitude a L");
        check(std::fabs(d2.tail_amplitude()-0.2)<1.0e-14, "type 2 tail amplitude a L");
        check(std::fabs(d3.tail_amplitude()-0.2)<1.0e-14, "type 3 tail amplitude a L");
        // carangiform: the head amplitude 0.02 L (a = 0.1), the minimum near s = 0.25
        d2.set_time(0.25);
        check(std::fabs(std::fabs(d2.disp(d2.xhead())) - 0.04)<1.0e-14, "type 2 head amplitude 0.02 L");
    }

    std::cout<<"deform_test: "<<npass<<" passed, "<<nfail<<" failed"<<std::endl;
    return nfail==0 ? 0 : 1;
}
