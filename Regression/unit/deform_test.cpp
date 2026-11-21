// Standalone check of the prescribed body deformation of a 6DOF body (X 195, src/6DOF_deform.cpp):
// no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY deform_test.cpp -o deform_test
// Run:    ./deform_test
// Checks: mass distribution along x from the closed hull (divergence theorem) for a box and for the
// 2D fish profile of Kern & Koumoutsakos (2006) extruded in y (volume, bin volumes against the
// width profile, also with the hull oriented inwards); head, tail and length; the envelopes
// (tail amplitudes of types 1-3); the deformation velocity is the time derivative of the
// displacement (with and without the ramp); recoil: no net lateral momentum and no angular momentum
// of the deformation (displacement and velocity) about the CoG, and without recoil the raw wave;
// inextensible midline (types 11-13): length of the midline, recoil (axial and lateral momentum,
// angular momentum), velocity = dd/dt, the velocity at the deformed point from the inverse map; the
// discrete recoil correction (X 196 2 with X 198, types 1 and 11): velocity = dd/dt and the inverse map with it;
// gait control (X 199): the state follows its targets, constant targets give X 195, the curvature
// bias, the bias along the envelope (X 200 w 1), and with changing targets midline length, recoil,
// velocity = dd/dt and the inverse map.
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

    // ---- inextensible midline (types 11-13): midline length, recoil, velocity = dd/dt, inverse map
    {
        std::vector<double> xs;
        const double b = 0.25;
        hull hf = fish2D(400,b,xs);
        double A = 0.0, Ax = 0.0;
        for(size_t i=0; i+1<xs.size(); ++i)
        {
            const double wa = fish_w(1.0-xs[i]), wb = fish_w(1.0-xs[i+1]);
            const double dx = xs[i+1]-xs[i];
            A  += dx*(wa + wb);
            Ax += dx*(wa + wb)*0.5*(xs[i]+xs[i+1]);
        }
        const double xc = Ax/A;
        
        for(int type : {11, 12, 13})
        {
            sixdof_geometry g(0);
            hf.to(g,xc,false);
            sixdof_deform d(type, type==11 ? 0.125 : 0.1, 1.0, 1.0, 1.0, 2, true);
            d.setup(g);
            check(d.inextensible(), "type "+std::to_string(type)+" inextensible");
            
            for(double t : {0.4, 1.3, 2.71})
            {
                d.set_time(t);
                
                // length of the deformed midline (l = 0)
                double len = 0.0, px = 0.0, pl = 0.0;
                const int K = 4000;
                for(int k=0; k<=K; ++k)
                {
                    const double x = d.xtail() + k*d.length()/K;
                    double dx, dl;
                    d.displacement(x,0.0,dx,dl);
                    if(k>0) len += std::hypot(x + dx - px, dl - pl);
                    px = x + dx; pl = dl;
                }
                check(std::fabs(len - 1.0)<2.0e-4, "type "+std::to_string(type)+" midline length "+std::to_string(len)+" at t = "+std::to_string(t));
                
                // momentum of the displacement and of the velocity, angular momentum of the velocity
                // about the CoG (the exact width; the recoil rotation is integrated in time)
                double Dx = 0.0, Dl = 0.0, Lz = 0.0, Vx = 0.0, Vl = 0.0, Lt = 0.0, href = 0.0, vref = 0.0;
                const int M = 4000, Q = 16;
                for(int k=0; k<M; ++k)
                {
                    const double xx = (k+0.5)/double(M);
                    const double w = fish_w(1.0-xx);
                    for(int q=0; q<Q; ++q)
                    {
                        const double l = -w + (q+0.5)*2.0*w/Q;
                        const double m = (2.0*w/Q)*b/M;
                        const double x = xx - xc;
                        double dx, dl, vx, vl;
                        d.displacement(x,l,dx,dl);
                        d.velocity(x,l,vx,vl);
                        // the angular momentum about the CoG with the deformed position
                        Dx += m*dx; Dl += m*dl; Lz += 0.0;
                        Vx += m*vx; Vl += m*vl; Lt += m*((x + dx)*vl - (l + dl)*vx);
                        href += m*std::hypot(dx,dl);
                        vref += m*std::hypot(vx,vl);
                    }
                }
                check(std::fabs(Dx)<3.0e-3*href && std::fabs(Dl)<3.0e-3*href && std::fabs(Lz)<3.0e-3*href,
                      "type "+std::to_string(type)+" recoil, displacement, t = "+std::to_string(t)+": "+std::to_string(Dx/href)+" "+std::to_string(Dl/href)+" "+std::to_string(Lz/href));
                check(std::fabs(Vx)<3.0e-3*vref && std::fabs(Vl)<3.0e-3*vref && std::fabs(Lt)<3.0e-3*vref,
                      "type "+std::to_string(type)+" recoil, velocity, t = "+std::to_string(t)+": "+std::to_string(Vx/vref)+" "+std::to_string(Vl/vref)+" "+std::to_string(Lt/vref));
            }
            
            // velocity = dd/dt of the material points, and the velocity at the deformed position
            for(double t : {0.3, 0.7, 1.6, 3.2})
            {
                double emax = 0.0, eimax = 0.0, vmax = 0.0;
                const double dt = 1.0e-6;
                for(int k=0; k<=40; ++k)
                for(double l : {-0.03, 0.0, 0.02})
                {
                    const double x = d.xtail() + 0.002 + k*(d.length()-0.004)/40.0;
                    if(std::fabs(l) > fish_w(1.0-(x+xc)) ) continue;
                    double dxp, dlp, dxm, dlm, vx, vl, ux, ul, dx, dl;
                    d.set_time(t+dt); d.displacement(x,l,dxp,dlp);
                    d.set_time(t-dt); d.displacement(x,l,dxm,dlm);
                    d.set_time(t);
                    d.velocity(x,l,vx,vl);
                    d.displacement(x,l,dx,dl);
                    d.velocity_at(x+dx,l+dl,ux,ul);
                    emax = std::max(emax, std::hypot((dxp-dxm)/(2.0*dt) - vx, (dlp-dlm)/(2.0*dt) - vl));
                    eimax = std::max(eimax, std::hypot(ux - vx, ul - vl));
                    vmax = std::max(vmax, std::hypot(vx,vl));
                }
                check(emax<1.0e-5*vmax, "type "+std::to_string(type)+" velocity = dd/dt at t = "+std::to_string(t)+": "+std::to_string(emax/vmax));
                check(eimax<1.0e-6*vmax, "type "+std::to_string(type)+" velocity at the deformed point at t = "+std::to_string(t)+": "+std::to_string(eimax/vmax));
            }
        }
        
        // without recoil: the lateral tail amplitude a L, the axial contraction of the tail
        sixdof_geometry g(0);
        hf.to(g,xc,false);
        sixdof_deform d(11, 0.125, 1.0, 1.0, 0.0, 2, false);
        d.setup(g);
        double hmax = 0.0, cmin = 1.0, cmax = 0.0;
        for(int k=0; k<1000; ++k)
        {
            d.set_time(k/1000.0);
            double dx, dl;
            d.displacement(d.xtail(),0.0,dx,dl);
            hmax = std::max(hmax, std::fabs(dl));
            cmin = std::min(cmin, dx);
            cmax = std::max(cmax, dx);
        }
        check(std::fabs(hmax-0.125)<2.0e-3, "type 11 tail amplitude: "+std::to_string(hmax));
        check(cmin>0.0 && cmax<0.1, "type 11 axial contraction at the tail in (0, 0.1): "+std::to_string(cmin)+" "+std::to_string(cmax));
    }
    
    // ---- discrete recoil correction (X 198): a rigid motion in the plane after the recoil
    {
        std::vector<double> xs;
        hull hf = fish2D(400,0.25,xs);
        for(int type : {1, 11})
        {
            sixdof_geometry g(0);
            hf.to(g,0.45,false);
            sixdof_deform d(type, 0.125, 1.0, 1.0, 1.0, 2, true);
            d.setup(g);
            
            // the correction moves linearly in time: tx = 0.01 + 0.02 t, tl = -0.005 + 0.03 t, psi = 0.02 - 0.05 t
            auto corr = [&](double t)
            {
                sixdof_deform::correction c;
                c.tx = 0.01 + 0.02*t;  c.vx = 0.02;
                c.tl = -0.005 + 0.03*t; c.vl = 0.03;
                c.psi = 0.02 - 0.05*t;  c.psit = -0.05;
                return c;
            };
            
            for(double t : {0.6, 1.7})
            {
                double emax = 0.0, eimax = 0.0, vmax = 0.0;
                const double dt = 1.0e-6;
                for(int k=0; k<=40; ++k)
                for(double l : {-0.02, 0.0, 0.015})
                {
                    const double x = d.xtail() + 0.003 + k*(d.length()-0.006)/40.0;
                    if(std::fabs(l) > fish_w(1.0-(x+0.45))) continue;
                    double dxp, dlp, dxm, dlm, vx, vl, ux, ul, dx, dl;
                    d.set_time(t+dt); d.set_corr(corr(t+dt)); d.displacement(x,l,dxp,dlp);
                    d.set_time(t-dt); d.set_corr(corr(t-dt)); d.displacement(x,l,dxm,dlm);
                    d.set_time(t); d.set_corr(corr(t));
                    d.velocity(x,l,vx,vl);
                    d.displacement(x,l,dx,dl);
                    d.velocity_at(x+dx,l+dl,ux,ul);
                    emax = std::max(emax, std::hypot((dxp-dxm)/(2.0*dt) - vx, (dlp-dlm)/(2.0*dt) - vl));
                    eimax = std::max(eimax, std::hypot(ux - vx, ul - vl));
                    vmax = std::max(vmax, std::hypot(vx,vl));
                }
                check(emax<1.0e-5*vmax, "type "+std::to_string(type)+" corrected velocity = dd/dt at t = "+std::to_string(t)+": "+std::to_string(emax/vmax));
                check(eimax<1.0e-6*vmax, "type "+std::to_string(type)+" corrected velocity at the deformed point at t = "+std::to_string(t)+": "+std::to_string(eimax/vmax));
            }
        }
    }
    
    // ---- gait control (X 199): the state follows its targets (95 % in t_s, continuous), the phase
    // is 2 pi int f; constant targets reproduce X 195; the curvature bias; and for the inextensible
    // midline with targets changing: midline length, recoil, velocity = dd/dt, the inverse map
    {
        std::vector<double> xs;
        const double b = 0.25;
        hull hf = fish2D(400,b,xs);
        double A = 0.0, Ax = 0.0;
        for(size_t i=0; i+1<xs.size(); ++i)
        {
            const double wa = fish_w(1.0-xs[i]), wb = fish_w(1.0-xs[i+1]);
            const double dx = xs[i+1]-xs[i];
            A  += dx*(wa + wb);
            Ax += dx*(wa + wb)*0.5*(xs[i]+xs[i+1]);
        }
        const double xc = Ax/A;
        
        const std::vector<sixdof_deform::gait_target> tg = {{0.5, 1.5, 0.8, 0.3, 0.2}, {1.0, 0.8, 1.2, -0.5, 0.0}};
        
        // the state: 95 % of the step from rest in t_s, continuity at the target times
        {
            sixdof_geometry g(0);
            hf.to(g,xc,false);
            sixdof_deform d(11, 0.125, 1.0, 1.0, 1.0, 2, true);
            d.setup(g);
            d.set_gait({{0.5, 1.5, 0.8, 0.3, 0.2}});
            d.set_time(0.7);
            check(std::fabs((d.frequency() - 1.5)/(1.0 - 1.5) - 0.05)<1.0e-3, "gait: f 95 % in t_s: "+std::to_string(d.frequency()));
            check(std::fabs((d.amp_factor() - 0.8)/(1.0 - 0.8) - 0.05)<1.0e-3, "gait: a 95 % in t_s: "+std::to_string(d.amp_factor()));
            check(std::fabs((d.bias() - 0.3)/(0.0 - 0.3) - 0.05)<1.0e-3, "gait: b 95 % in t_s: "+std::to_string(d.bias()));
            d.set_time(0.3);
            check(d.frequency()==1.0 && d.amp_factor()==1.0 && d.bias()==0.0, "gait: start values before the first target");
            
            d.set_gait(tg);
            double fm, fp, bm, bp;
            d.set_time(1.0 - 1.0e-9); fm = d.frequency(); bm = d.bias();
            d.set_time(1.0 + 1.0e-9); fp = d.frequency(); bp = d.bias();
            check(std::fabs(fp - fm)<1.0e-7 && std::fabs(bp - bm)<1.0e-7, "gait: continuous at a target time");
            d.set_time(6.0);
            check(std::fabs(d.frequency() - 0.8)<1.0e-9 && std::fabs(d.amp_factor() - 1.2)<1.0e-9 && std::fabs(d.bias() + 0.5)<1.0e-9, "gait: settled on the last targets");
        }
        
        // constant targets equal to the X 195 values: the same body (types 1 and 11)
        for(int type : {1, 11})
        {
            sixdof_geometry g(0), g2(0);
            hf.to(g,xc,false);
            hf.to(g2,xc,false);
            sixdof_deform d(type, 0.125, 1.0, 1.0, 1.0, 2, true), e(type, 0.125, 1.0, 1.0, 1.0, 2, true);
            d.setup(g);
            e.setup(g2);
            e.set_gait({{2.0, 1.0, 1.0, 0.0, 0.0}});
            double emax = 0.0;
            for(double t : {0.4, 1.3, 2.5})
            {
                d.set_time(t); e.set_time(t);
                for(int k=0; k<=40; ++k)
                {
                    const double x = d.xtail() + k*d.length()/40.0;
                    double dx, dl, ex, el;
                    d.displacement(x,0.01,dx,dl);
                    e.displacement(x,0.01,ex,el);
                    emax = std::max(emax, std::hypot(dx-ex, dl-el));
                }
            }
            check(emax<1.0e-12, "type "+std::to_string(type)+" gait with the X 195 values = X 195: "+std::to_string(emax));
        }
        
        // the curvature bias alone (a = 0, no recoil, type 1): h = b L s^2/2, curvature b/L
        {
            sixdof_geometry g(0);
            hf.to(g,xc,false);
            sixdof_deform d(1, 0.125, 1.0, 1.0, 0.0, 2, false);
            d.setup(g);
            d.set_gait({{0.0, 1.0, 0.0, 0.4, 0.0}});
            d.set_time(0.37);
            double emax = 0.0;
            for(double s : {0.2, 0.5, 0.9})
            emax = std::max(emax, std::fabs(d.disp(d.xhead() - s) - 0.4*s*s/2.0));
            check(emax<1.0e-14, "curvature bias b L s^2/2: "+std::to_string(emax));
        }
        
        // the bias along the envelope (X 200 w 1; a = 0, no recoil, type 12): the angle of the midline
        // b int_0^s A / int_0^1 A, the curvature with the carangiform envelope (s from the head)
        {
            sixdof_geometry g(0);
            hf.to(g,xc,false);
            sixdof_deform d(12, 0.075, 0.75, 1.0, 0.0, 2, false);
            d.setup(g);
            d.set_bias_shape(1);
            d.set_gait({{0.0, 1.0, 0.0, 0.4, 0.0}});
            d.set_time(0.37);
            auto I = [](double x) {return 0.02*x - 0.04125*x*x + 0.1625/3.0*x*x*x;};
            double emax = 0.0;
            const double h = 1.0e-4;
            for(double sg : {0.1, 0.3, 0.5, 0.7, 0.9, 0.99})
            {
                double dx1, dl1, dx2, dl2;
                const double x1 = d.xhead() - (sg - h)*d.length(), x2 = d.xhead() - (sg + h)*d.length();
                d.displacement(x1,0.0,dx1,dl1);
                d.displacement(x2,0.0,dx2,dl2);
                const double th = std::atan2(dl2 - dl1, (x1 + dx1) - (x2 + dx2));
                emax = std::max(emax, std::fabs(th - 0.4*I(sg)/I(1.0)));
            }
            check(emax<5.0e-4, "type 12, the bias along the envelope: the angle b int A / int A, error "+std::to_string(emax));
        }
        
        // inextensible midline with changing targets
        for(int type : {1, 11, 12})
        {
            sixdof_geometry g(0);
            hf.to(g,xc,false);
            sixdof_deform d(type, 0.125, 1.0, 1.0, 0.5, 2, true);
            d.setup(g);
            if(type==12)
            d.set_bias_shape(1);
            d.set_gait(tg);
            const std::string ty = "type "+std::to_string(type)+" gait"+(type==12 ? " (the bias along the envelope)" : "");
            
            for(double t : {0.6, 0.9, 1.1, 1.4, 2.2})
            {
                d.set_time(t);
                
                if(type==11 || type==12)
                {
                    double len = 0.0, px = 0.0, pl = 0.0;
                    const int K = 4000;
                    for(int k=0; k<=K; ++k)
                    {
                        const double x = d.xtail() + k*d.length()/K;
                        double dx, dl;
                        d.displacement(x,0.0,dx,dl);
                        if(k>0) len += std::hypot(x + dx - px, dl - pl);
                        px = x + dx; pl = dl;
                    }
                    check(std::fabs(len - 1.0)<2.0e-4, ty+" midline length "+std::to_string(len)+" at t = "+std::to_string(t));
                    
                    double Dx = 0.0, Dl = 0.0, Vx = 0.0, Vl = 0.0, Lt = 0.0, href = 0.0, vref = 0.0;
                    const int M = 4000, Q = 16;
                    for(int k=0; k<M; ++k)
                    {
                        const double xx = (k+0.5)/double(M);
                        const double w = fish_w(1.0-xx);
                        for(int q=0; q<Q; ++q)
                        {
                            const double l = -w + (q+0.5)*2.0*w/Q;
                            const double m = (2.0*w/Q)*b/M;
                            const double x = xx - xc;
                            double dx, dl, vx, vl;
                            d.displacement(x,l,dx,dl);
                            d.velocity(x,l,vx,vl);
                            Dx += m*dx; Dl += m*dl;
                            Vx += m*vx; Vl += m*vl; Lt += m*((x + dx)*vl - (l + dl)*vx);
                            href += m*std::hypot(dx,dl);
                            vref += m*std::hypot(vx,vl);
                        }
                    }
                    check(std::fabs(Dx)<3.0e-3*href && std::fabs(Dl)<3.0e-3*href, ty+" recoil, displacement, t = "+std::to_string(t)+": "+std::to_string(Dx/href)+" "+std::to_string(Dl/href));
                    check(std::fabs(Vx)<3.0e-3*vref && std::fabs(Vl)<3.0e-3*vref && std::fabs(Lt)<3.0e-3*vref,
                          ty+" recoil, velocity, t = "+std::to_string(t)+": "+std::to_string(Vx/vref)+" "+std::to_string(Vl/vref)+" "+std::to_string(Lt/vref));
                }
                
                // velocity = dd/dt and the velocity at the deformed point
                double emax = 0.0, eimax = 0.0, vmax = 0.0;
                const double dt = 1.0e-6;
                for(int k=0; k<=40; ++k)
                for(double l : {-0.03, 0.0, 0.02})
                {
                    const double x = d.xtail() + 0.002 + k*(d.length()-0.004)/40.0;
                    if(std::fabs(l) > fish_w(1.0-(x+xc)) ) continue;
                    double dxp, dlp, dxm, dlm, vx, vl, ux, ul, dx, dl;
                    d.set_time(t+dt); d.displacement(x,l,dxp,dlp);
                    d.set_time(t-dt); d.displacement(x,l,dxm,dlm);
                    d.set_time(t);
                    d.velocity(x,l,vx,vl);
                    d.displacement(x,l,dx,dl);
                    d.velocity_at(x+dx,l+dl,ux,ul);
                    emax = std::max(emax, std::hypot((dxp-dxm)/(2.0*dt) - vx, (dlp-dlm)/(2.0*dt) - vl));
                    eimax = std::max(eimax, std::hypot(ux - vx, ul - vl));
                    vmax = std::max(vmax, std::hypot(vx,vl));
                }
                check(emax<1.0e-5*vmax, ty+" velocity = dd/dt at t = "+std::to_string(t)+": "+std::to_string(emax/vmax));
                check(eimax<1.0e-6*vmax, ty+" velocity at the deformed point at t = "+std::to_string(t)+": "+std::to_string(eimax/vmax));
            }
        }
    }
    
    // ---- gait control: a C-start (bias 3, the tail at 172 deg to the head) on the inextensible
    // midline: length, recoil, velocity = dd/dt and the inverse map
    {
        std::vector<double> xs;
        const double b = 0.25;
        hull hf = fish2D(400,b,xs);
        double A = 0.0, Ax = 0.0;
        for(size_t i=0; i+1<xs.size(); ++i)
        {
            const double wa = fish_w(1.0-xs[i]), wb = fish_w(1.0-xs[i+1]);
            const double dx = xs[i+1]-xs[i];
            A  += dx*(wa + wb);
            Ax += dx*(wa + wb)*0.5*(xs[i]+xs[i+1]);
        }
        const double xc = Ax/A;
        sixdof_geometry g(0);
        hf.to(g,xc,false);
        sixdof_deform d(11, 0.125, 1.0, 1.0, 0.0, 2, true);
        d.setup(g);
        d.set_gait({{0.2, 1.0, 0.3, 3.0, 0.1}});
        
        for(double t : {0.25, 0.32, 0.6})
        {
            d.set_time(t);
            double len = 0.0, px = 0.0, pl = 0.0;
            const int K = 4000;
            for(int k=0; k<=K; ++k)
            {
                const double x = d.xtail() + k*d.length()/K;
                double dx, dl;
                d.displacement(x,0.0,dx,dl);
                if(k>0) len += std::hypot(x + dx - px, dl - pl);
                px = x + dx; pl = dl;
            }
            check(std::fabs(len - 1.0)<2.0e-4, "C-start midline length "+std::to_string(len)+" at t = "+std::to_string(t));
            
            double Vx = 0.0, Vl = 0.0, Lt = 0.0, vref = 0.0, emax = 0.0, eimax = 0.0, vmax = 0.0;
            const int M = 2000, Q = 8;
            for(int k=0; k<M; ++k)
            {
                const double xx = (k+0.5)/double(M);
                const double w = fish_w(1.0-xx);
                for(int q=0; q<Q; ++q)
                {
                    const double l = -w + (q+0.5)*2.0*w/Q;
                    const double m = (2.0*w/Q)*b/M;
                    const double x = xx - xc;
                    double dx, dl, vx, vl;
                    d.displacement(x,l,dx,dl);
                    d.velocity(x,l,vx,vl);
                    Vx += m*vx; Vl += m*vl; Lt += m*((x + dx)*vl - (l + dl)*vx);
                    vref += m*std::hypot(vx,vl);
                }
            }
            check(std::fabs(Vx)<3.0e-3*vref && std::fabs(Vl)<3.0e-3*vref && std::fabs(Lt)<3.0e-3*vref,
                  "C-start recoil, velocity, t = "+std::to_string(t)+": "+std::to_string(Vx/vref)+" "+std::to_string(Vl/vref)+" "+std::to_string(Lt/vref));
            
            const double dt = 1.0e-6;
            for(int k=0; k<=40; ++k)
            for(double l : {-0.03, 0.0, 0.02})
            {
                const double x = d.xtail() + 0.002 + k*(d.length()-0.004)/40.0;
                if(std::fabs(l) > fish_w(1.0-(x+xc)) ) continue;
                double dxp, dlp, dxm, dlm, vx, vl, ux, ul, dx, dl;
                d.set_time(t+dt); d.displacement(x,l,dxp,dlp);
                d.set_time(t-dt); d.displacement(x,l,dxm,dlm);
                d.set_time(t);
                d.velocity(x,l,vx,vl);
                d.displacement(x,l,dx,dl);
                d.velocity_at(x+dx,l+dl,ux,ul);
                emax = std::max(emax, std::hypot((dxp-dxm)/(2.0*dt) - vx, (dlp-dlm)/(2.0*dt) - vl));
                eimax = std::max(eimax, std::hypot(ux - vx, ul - vl));
                vmax = std::max(vmax, std::hypot(vx,vl));
            }
            check(emax<1.0e-5*vmax, "C-start velocity = dd/dt at t = "+std::to_string(t)+": "+std::to_string(emax/vmax));
            check(eimax<1.0e-6*vmax, "C-start velocity at the deformed point at t = "+std::to_string(t)+": "+std::to_string(eimax/vmax));
        }
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
