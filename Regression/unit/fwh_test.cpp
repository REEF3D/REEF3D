// Standalone verification of the permeable FW-H kernel (fwh_permeable): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src fwh_test.cpp ../../src/acoustics_fwh_kernel.cpp -o fwh_test
// Run:    ./fwh_test
//
// The surface data are the exact fields of acoustic point monopoles (linear acoustics), so the
// FW-H integral has to reproduce the exact pressure outside the surface and zero inside.
// Tests: monopole and dipole pair at near and far observers in several directions, independence
// of the surface size, second-order convergence in panel size and time step, adaptive time steps,
// observer inside the surface, a source below a pressure-release surface (image observer, weight
// -1) against the image-source solution, panels split over two "ranks" summed as with MPI.
#include"acoustics_fwh_kernel.h"
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<random>
#include<cstdio>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static const double c0   = 1500.0;
static const double rho0 = 1000.0;
static const double A    = 1.0e-3;                 // volume flow amplitude [m^3/s]
static const double w1   = 2.0*M_PI*400.0;
static const double w2   = 2.0*M_PI*1500.0;        // shortest wavelength 1 m

// point monopole of volume flow s*q(t): phi = -s q(t-R/c0)/(4 pi R)
//  p' = rho0 s q'(tau)/(4 pi R),  u = s Rhat (q'(tau)/(4 pi R c0) + q(tau)/(4 pi R^2))
struct source
{
    double x[3];
    double s;
};

static double q(double t)  {return A*(sin(w1*t) + 0.5*sin(w2*t+0.7));}
static double qd(double t) {return A*(w1*cos(w1*t) + 0.5*w2*cos(w2*t+0.7));}

static void field(const std::vector<source> &src, const double *x, double t, double &p, double *u)
{
    p = u[0] = u[1] = u[2] = 0.0;

    for(const source &m : src)
    {
        const double r[3] = {x[0]-m.x[0], x[1]-m.x[1], x[2]-m.x[2]};
        const double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        const double tau = t - R/c0;

        p += rho0*m.s*qd(tau)/(4.0*M_PI*R);

        const double ur = m.s*(qd(tau)/(4.0*M_PI*R*c0) + q(tau)/(4.0*M_PI*R*R));
        for(int c=0; c<3; ++c)
        u[c] += ur*r[c]/R;
    }
}

// box [lo,hi] with n panels per metre (at least 1 per edge), midpoint panels, outward normals
static std::vector<fwh_panel> box(const double *lo, const double *hi, double n)
{
    std::vector<fwh_panel> pan;

    for(int a=0; a<3; ++a)
    for(int side=0; side<2; ++side)
    {
        const int b = (a+1)%3, c = (a+2)%3;
        const int nb = std::max(1,int(lround(n*(hi[b]-lo[b]))));
        const int nc = std::max(1,int(lround(n*(hi[c]-lo[c]))));
        const double db = (hi[b]-lo[b])/nb, dc = (hi[c]-lo[c])/nc;

        for(int i=0; i<nb; ++i)
        for(int j=0; j<nc; ++j)
        {
            fwh_panel P;
            P.x[a] = side ? hi[a] : lo[a];
            P.x[b] = lo[b] + (i+0.5)*db;
            P.x[c] = lo[c] + (j+0.5)*dc;
            P.n[0]=P.n[1]=P.n[2]=0.0;
            P.n[a] = side ? 1.0 : -1.0;
            P.dS = db*dc;
            pan.push_back(P);
        }
    }
    return pan;
}

struct obs_point
{
    double x[3];
    double img[3];      // image point (weight -1), used if mirror
};

struct setup
{
    std::vector<source> src;        // sources whose field is on the surface
    std::vector<source> exact;      // sources of the exact solution (with image sources)
    double lo[3], hi[3];
    double npm;                     // panels per metre
    double dt, dto, tend;
    double dtvar;                   // relative random variation of the time step (adaptive)
    bool mirror;
    int nsplit;                     // panels distributed round robin over nsplit "ranks"
    std::vector<obs_point> ob;
};

struct result
{
    std::vector<double> err;        // max|p - p_exact| / max|p_exact| over the complete range
    std::vector<double> amp;        // max|p| (for the observer inside)
    std::vector<std::vector<double>> sig;
    long k0;
};

static result run(const setup &S)
{
    std::vector<fwh_panel> all = box(S.lo,S.hi,S.npm);

    std::vector<fwh_permeable> rk(S.nsplit, fwh_permeable(c0,rho0,S.dto));
    std::vector<std::vector<int>> idx(S.nsplit);

    for(int i=0; i<int(all.size()); ++i)
    {
        rk[i%S.nsplit].add_panel(all[i]);
        idx[i%S.nsplit].push_back(i);
    }

    for(fwh_permeable &f : rk)
    for(const obs_point &o : S.ob)
    {
        const int id = f.add_observer(o.x);
        if(S.mirror)
        f.add_image(id,o.img,-1.0);
    }

    // mt19937 output is fixed by the standard (the distributions are not): same steps on all platforms
    std::mt19937 gen(12345);

    std::vector<double> p, u;
    double t = 0.0;
    while(t<=S.tend)
    {
        for(int r=0; r<S.nsplit; ++r)
        {
            const int np = int(idx[r].size());
            p.resize(np);
            u.resize(3*np);

            for(int i=0; i<np; ++i)
            field(S.src,all[idx[r][i]].x,t,p[i],&u[3*i]);

            rk[r].step(t,p.data(),u.data());
        }
        t += S.dt*(1.0 + S.dtvar*(2.0*double(gen())/double(gen.max()) - 1.0));
    }

    result R;
    R.k0 = rk[0].first_index();

    for(int o=0; o<int(S.ob.size()); ++o)
    {
        // sum over the ranks and global delay bounds, as the MPI coupling does
        std::vector<double> sig;
        double dmin=1.0e300, dmax=-1.0e300;

        for(fwh_permeable &f : rk)
        {
            const std::vector<double> &d = f.data(o);
            if(d.size()>sig.size())
            sig.resize(d.size(),0.0);
            for(size_t m=0; m<d.size(); ++m)
            sig[m] += d[m];

            double a, b;
            f.delays(o,a,b);
            dmin = std::min(dmin,a);
            dmax = std::max(dmax,b);
        }

        double t0, t1;
        if(!rk[0].complete_range(dmin,dmax,t0,t1))
        {
            R.err.push_back(1.0e30);
            R.amp.push_back(1.0e30);
            R.sig.push_back(sig);
            continue;
        }

        double emax=0.0, pmax=0.0, amax=0.0;
        for(size_t m=0; m<sig.size(); ++m)
        {
            const double tk = double(R.k0+long(m))*S.dto;
            if(tk<t0 || tk>t1)
            continue;

            double pe, ue[3];
            field(S.exact,S.ob[o].x,tk,pe,ue);

            emax = std::max(emax,fabs(sig[m]-pe));
            pmax = std::max(pmax,fabs(pe));
            amax = std::max(amax,fabs(sig[m]));
        }
        R.err.push_back(pmax>0.0 ? emax/pmax : emax);
        R.amp.push_back(amax);
        R.sig.push_back(sig);
    }
    return R;
}

static obs_point point(double x, double y, double z)
{
    return {{x,y,z},{x,y,-z}};
}

static setup base()
{
    setup S;
    S.src   = {{{0.0,0.0,0.0},1.0}};
    S.exact = S.src;
    for(int c=0; c<3; ++c)
    {
        S.lo[c] = -0.5;
        S.hi[c] =  0.5;
    }
    S.npm   = 40.0;
    S.dt    = 1.0/(1500.0*40.0);
    S.dto   = S.dt;
    S.tend  = 0.04;
    S.dtvar = 0.0;
    S.mirror = false;
    S.nsplit = 1;
    return S;
}

// peak pressure of the exact solution at x (for the observer inside)
static double exact_peak(const setup &S, const double *x)
{
    double pm = 0.0;
    for(double t=0.0; t<0.01; t+=S.dt)
    {
        double p, u[3];
        field(S.exact,x,t,p,u);
        pm = std::max(pm,fabs(p));
    }
    return pm;
}

int main()
{
    std::cout<<"FW-H permeable surface kernel"<<std::endl;

    // --- monopole, observers near and far in several directions
    {
        setup S = base();
        S.ob = {point(30.0,0.0,0.0), point(0.0,-30.0,0.0), point(17.0,17.0,17.0),
                point(1.5,0.0,0.0), point(1.2,1.2,-1.2)};
        result R = run(S);

        std::cout<<"monopole, 40 panels/m, 40 steps per period of 1500 Hz"<<std::endl;
        for(size_t o=0; o<S.ob.size(); ++o)
        {
            char s[200];
            snprintf(s,sizeof(s),"observer (%5.1f,%5.1f,%5.1f): relative error %.2e < 1e-2",
                     S.ob[o].x[0],S.ob[o].x[1],S.ob[o].x[2],R.err[o]);
            check(R.err[o]<1.0e-2,s);
        }
    }

    // --- dipole pair (two monopoles of opposite sign), directivity zero on the mid-plane
    {
        setup S = base();
        S.src   = {{{0.0,0.0,0.05},1.0},{{0.0,0.0,-0.05},-1.0}};
        S.exact = S.src;
        S.ob = {point(0.0,0.0,30.0), point(20.0,0.0,20.0), point(0.0,0.0,-1.5)};
        result R = run(S);

        std::cout<<"dipole pair"<<std::endl;
        for(size_t o=0; o<S.ob.size(); ++o)
        {
            char s[200];
            snprintf(s,sizeof(s),"observer (%5.1f,%5.1f,%5.1f): relative error %.2e < 1e-2",
                     S.ob[o].x[0],S.ob[o].x[1],S.ob[o].x[2],R.err[o]);
            check(R.err[o]<1.0e-2,s);
        }
    }

    // --- surface independence: two box sizes, off-centre source
    {
        std::cout<<"surface independence"<<std::endl;
        double e[2];
        for(int b=0; b<2; ++b)
        {
            setup S = base();
            S.src   = {{{0.1,-0.05,0.08},1.0}};
            S.exact = S.src;
            const double h = b==0 ? 0.4 : 1.0;
            for(int c=0; c<3; ++c)
            {
                S.lo[c] = -h;
                S.hi[c] =  h;
            }
            S.npm = b==0 ? 60.0 : 30.0;
            S.ob = {point(25.0,5.0,-3.0)};
            e[b] = run(S).err[0];
        }
        char s[200];
        snprintf(s,sizeof(s),"half width 0.4 m: %.2e, 1.0 m: %.2e, both < 1e-2",e[0],e[1]);
        check(e[0]<1.0e-2 && e[1]<1.0e-2,s);
    }

    // --- convergence in panel size (small time step) and time step (fine panels)
    {
        std::cout<<"convergence"<<std::endl;
        const double npm[3] = {5.0,10.0,20.0};
        double es[3];
        for(int l=0; l<3; ++l)
        {
            setup S = base();
            S.npm  = npm[l];
            S.dt   = 1.0/(1500.0*320.0);
            S.dto  = S.dt;
            S.tend = 0.03;
            S.ob = {point(1.5,0.3,0.2)};
            es[l] = run(S).err[0];
        }
        const double os = log2(es[1]/es[2]);
        char s[200];
        snprintf(s,sizeof(s),"panels 5/10/20 per m: %.2e %.2e %.2e, order %.2f in [1.8,2.3]",es[0],es[1],es[2],os);
        check(os>1.8 && os<2.3,s);

        const double spp[3] = {10.0,20.0,40.0};
        double et[3];
        for(int l=0; l<3; ++l)
        {
            setup S = base();
            S.npm  = 80.0;
            S.dt   = 1.0/(1500.0*spp[l]);
            S.dto  = 1.0/(1500.0*40.0);
            S.tend = 0.03;
            S.ob = {point(20.0,3.0,-2.0)};
            et[l] = run(S).err[0];
        }
        const double ot = log2(et[0]/et[1]);
        snprintf(s,sizeof(s),"steps 10/20/40 per period: %.2e %.2e %.2e, order %.2f in [1.8,2.3]",et[0],et[1],et[2],ot);
        check(ot>1.8 && ot<2.3,s);
    }

    // --- adaptive time steps (+-40 % random) as with the CFL control of REEF3D
    {
        setup S = base();
        S.dtvar = 0.4;
        S.ob = {point(30.0,0.0,0.0), point(1.5,0.0,0.0)};
        result R = run(S);
        char s[200];
        snprintf(s,sizeof(s),"adaptive time step: errors %.2e %.2e < 1.5e-2",R.err[0],R.err[1]);
        check(R.err[0]<1.5e-2 && R.err[1]<1.5e-2,s);
    }

    // --- observer inside the surface: the FW-H integral vanishes
    {
        setup S = base();
        S.ob = {point(0.2,-0.15,0.1)};
        result R = run(S);
        const double pe = exact_peak(S,S.ob[0].x);
        char s[200];
        snprintf(s,sizeof(s),"observer inside: max|p| / exact peak there = %.2e < 1e-2",R.amp[0]/pe);
        check(R.amp[0]/pe<1.0e-2,s);
    }

    // --- pressure-release free surface z=0: source at depth 1 m, box below the surface,
    //     image observers with weight -1 against the exact image-source solution
    {
        setup S = base();
        S.src = {{{0.0,0.0,-1.0},1.0}};
        S.exact = {{{0.0,0.0,-1.0},1.0},{{0.0,0.0,1.0},-1.0}};
        S.lo[2] = -1.5;
        S.hi[2] = -0.5;
        S.mirror = true;
        S.ob = {point(25.0,0.0,-5.0), point(3.0,2.0,-0.3), point(0.0,0.0,-30.0)};
        result R = run(S);

        std::cout<<"free surface mirror (Lloyd's mirror)"<<std::endl;
        for(size_t o=0; o<S.ob.size(); ++o)
        {
            char s[200];
            snprintf(s,sizeof(s),"observer (%5.1f,%5.1f,%5.1f): relative error %.2e < 1e-2",
                     S.ob[o].x[0],S.ob[o].x[1],S.ob[o].x[2],R.err[o]);
            check(R.err[o]<1.0e-2,s);
        }
    }

    // --- panels split over 3 "ranks", summed: same signal as one rank
    {
        setup S = base();
        S.npm  = 20.0;
        S.tend = 0.03;
        S.ob = {point(10.0,4.0,-2.0)};
        result R1 = run(S);
        S.nsplit = 3;
        result R3 = run(S);

        double dmax=0.0, pmax=0.0;
        const size_t n = std::min(R1.sig[0].size(),R3.sig[0].size());
        for(size_t m=0; m<n; ++m)
        {
            dmax = std::max(dmax,fabs(R1.sig[0][m]-R3.sig[0][m]));
            pmax = std::max(pmax,fabs(R1.sig[0][m]));
        }
        char s[200];
        snprintf(s,sizeof(s),"3 ranks summed vs 1 rank: same first index, relative difference %.2e < 1e-12",dmax/pmax);
        check(R1.k0==R3.k0 && R1.sig[0].size()==R3.sig[0].size() && dmax/pmax<1.0e-12,s);
    }

    std::cout<<(nfail ? "FAILED: " : "all passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
