// Standalone verification of the REEF3D rod-tree solver (src/rodtree.cpp).
// Build: make -C tests/rodtree && ./tests/rodtree/test_rodtree
//
//  1  cantilever under self-weight vs Euler-Bernoulli, N-convergence
//  2  first bending frequency vs analytical (explicit and implicit)
//  3  objectivity: rigid rotation of a branching colony is stress free
//  4  energy conservation of the explicit integrator (branching, undamped)
//  5  drag reconfiguration: implicit vs explicit steady state
//  6  branching colony in oscillatory flow: implicit vs explicit
//  7  large-deflection cantilever vs elastica (tip load via drag-free end mass)

#include"../../src/rodtree.h"
#include<cstdio>
#include<cmath>
#include<sstream>
#include<vector>
#include<string>

typedef Eigen::Vector3d V3;
static const double PI = 3.14159265358979323846;
static int failures = 0;

static void check(bool ok, const char* what)
{
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what);
    if(!ok) ++failures;
}

static std::string cantilever(double L, double r, double E, double rho, int nref, V3 dir, double beta, double Cdn=1.2, double Ca=1.0)
{
    std::ostringstream s;
    s<<"colony c\nmaterial "<<E<<" 0.4 "<<rho<<"\nhydro "<<Cdn<<" 0.0 "<<Ca<<"\ndamping "<<beta<<"\nrefine "<<nref<<"\n";
    s<<"node 1 0 0 0 "<<r<<"\nnode 2 "<<L*dir(0)<<" "<<L*dir(1)<<" "<<L*dir(2)<<" "<<r<<"\nedge 1 2\nclamp 1\nend\n";
    return s.str();
}

// simple gorgonian-like fan: stem + 2 levels of branches
static std::string colony_fan(double beta, int nref)
{
    std::ostringstream s;
    s<<"colony fan\nmaterial 5.0e6 0.4 1100\nhydro 1.2 0.01 1.0\ndamping "<<beta<<"\nrefine "<<nref<<"\n";
    s<<"node 1 0 0 0 0.006\nnode 2 0 0 0.10 0.005\n"
     <<"node 3 0.05 0 0.18 0.004\nnode 4 -0.05 0 0.18 0.004\nnode 5 0 0.02 0.20 0.004\n"
     <<"node 6 0.08 0 0.26 0.003\nnode 7 0.03 0.01 0.27 0.003\nnode 8 -0.08 0 0.26 0.003\nnode 9 -0.03 -0.01 0.27 0.003\n";
    s<<"edge 1 2\nedge 2 3\nedge 2 4\nedge 2 5\nedge 3 6\nedge 3 7\nedge 4 8\nedge 4 9\nclamp 1\nend\n";
    return s.str();
}

static void load(rodtree& rt, const std::string& in)
{
    std::istringstream is(in);
    rt.read(is);
    rt.finalize_setup();
}

static void set_fluid(rodtree& rt, const V3& u, const V3& a, double chi)
{
    for(int e=0; e<rt.nelem(); ++e)
    {
        rt.elem(e).uf = u;
        rt.elem(e).af = a;
        rt.elem(e).chi = chi;
    }
}

int main()
{
    setvbuf(stdout,nullptr,_IONBF,0);
    // common beam
    const double L = 0.3, r = 0.004, E = 5.0e6, rho = 1100.0;
    const double A = PI*r*r, I = 0.25*PI*r*r*r*r, EI = E*I;

    // ------------------------------------------------------------------ 1
    std::printf("Test 1: cantilever under self-weight (small deflection)\n");
    {
        const double gz = -0.05;                    // small load -> linear regime
        const double q = rho*A*std::fabs(gz);
        const double w_ref = q*L*L*L*L/(8.0*EI);
        double err_prev = 0.0;
        for(int n : {5,10,20,40})
        {
            rodtree rt;
            load(rt,cantilever(L,r,E,rho,n,V3(1,0,0),0.4));
            rt.set_gravity(V3(0,0,gz));
            rt.set_fluid_density(0.0);
            for(int s=0; s<5000; ++s) rt.advance(1.0e-3);
            double w = -rt.tip_displacement(0)(2);
            double err = std::fabs(w-w_ref)/w_ref;
            std::printf("    N=%3d  tip %.6e  ref %.6e  rel.err %.3e", n, w, w_ref, err);
            if(err_prev>0.0) std::printf("  ratio %.2f", err_prev/err);
            std::printf("\n");
            if(n==40) check(err<2.0e-3,"N=40 within 0.2% of q L^4/(8EI)");
            err_prev = err;
        }
    }

    // ------------------------------------------------------------------ 2
    std::printf("Test 2: first bending frequency, vertical cantilever in vacuum\n");
    {
        const double f_ref = (1.87510407*1.87510407/(2.0*PI))*std::sqrt(EI/(rho*A*L*L*L*L));
        for(int integ=0; integ<2; ++integ)
        {
            rodtree rt;
            load(rt,cantilever(L,r,E,rho,40,V3(0,0,1),0.0));
            rt.set_fluid_density(0.0);
            rt.set_integrator(0);                   // damped pre-deflection: implicit
            // static pre-deflection by lateral gravity (damped phase)
            rt.set_gravity(V3(0.05,0,0));
            for(int e=0; e<rt.njoint(); ++e) const_cast<rodtree::joint&>(rt.jnt(e)).beta = 0.4;
            for(int s=0; s<3000; ++s) rt.advance(1.0e-3);
            for(int e=0; e<rt.njoint(); ++e) const_cast<rodtree::joint&>(rt.jnt(e)).beta = 0.0;
            for(int e=0; e<rt.nelem(); ++e) {rt.elem(e).v.setZero(); rt.elem(e).w.setZero();}
            rt.set_gravity(V3(0,0,0));
            rt.set_integrator(integ);
            const double dt = integ==0 ? 2.0e-4 : 1.0e-3;
            if(integ==0) rt.set_substeps(1);
            double prev = rt.tip_displacement(0)(0), tprev = 0.0, t = 0.0;
            std::vector<double> cross;
            double amp0 = std::fabs(prev), ampN = 0.0;
            while(t<6.0/f_ref)
            {
                rt.advance(dt); t += dt;
                double x = rt.tip_displacement(0)(0);
                if(prev<0.0 && x>=0.0) cross.push_back(tprev + dt*(-prev)/(x-prev));
                ampN = std::max(ampN,std::fabs(x)*(t>5.0/f_ref));
                prev = x; tprev = t;
            }
            double f = (cross.size()>1) ? (cross.size()-1)/(cross.back()-cross.front()) : 0.0;
            double err = std::fabs(f-f_ref)/f_ref;
            std::printf("    %s: f %.5f Hz  ref %.5f Hz  rel.err %.3e  amplitude ratio after 5 periods %.4f  (substeps/step %d)\n",
                        integ==0 ? "implicit" : "explicit", f, f_ref, err, ampN/amp0, rt.substeps_used());
            check(err<5.0e-3, integ==0 ? "implicit frequency within 0.5%" : "explicit frequency within 0.5%");
        }
    }

    // ------------------------------------------------------------------ 3
    std::printf("Test 3: objectivity of a branching colony under rigid rotation\n");
    {
        rodtree rt;
        load(rt,colony_fan(0.0,4));
        Eigen::Quaterniond R(Eigen::AngleAxisd(0.9,V3(0.3,-0.5,0.8).normalized()));
        for(int e=0; e<rt.nelem(); ++e)
        {
            rt.elem(e).c = R*rt.elem(e).c;
            rt.elem(e).q = R*rt.elem(e).q;
        }
        // anchors live in the world, so rotate them as well
        for(int j=0; j<rt.njoint(); ++j)
        if(rt.jnt(j).a<0) const_cast<rodtree::joint&>(rt.jnt(j)).pa = R*rt.jnt(j).pa;
        // the clamp also fixes orientation: rotate its rest frame
        for(int j=0; j<rt.njoint(); ++j)
        if(rt.jnt(j).a<0) const_cast<rodtree::joint&>(rt.jnt(j)).qrel0 = R*rt.jnt(j).qrel0;
        double Ee = rt.elastic_energy();
        std::printf("    elastic energy after rigid rotation: %.3e J\n", Ee);
        check(Ee<1.0e-18,"rigidly rotated colony is stress free");
    }

    // ------------------------------------------------------------------ 4
    std::printf("Test 4: explicit integrator energy conservation (undamped fan)\n");
    {
        rodtree rt;
        load(rt,colony_fan(0.0,4));
        rt.set_gravity(V3(0,0,0));
        rt.set_fluid_density(0.0);
        rt.set_integrator(1);
        for(int e=0; e<rt.nelem(); ++e)
        rt.elem(e).v = V3(0.3,0.1,0.0)*(rt.elem(e).c(2)/0.27);
        double E0 = rt.kinetic_energy() + rt.elastic_energy(), Emax = 0.0;
        for(int s=0; s<1000; ++s)
        {
            rt.advance(1.0e-3);
            Emax = std::max(Emax,std::fabs(rt.kinetic_energy()+rt.elastic_energy()-E0)/E0);
        }
        std::printf("    E0 %.4e J  max |dE|/E0 over 1 s: %.3e  (substeps/step %d)\n", E0, Emax, rt.substeps_used());
        check(Emax<2.0e-2,"energy drift below 2% (symplectic, bounded)");
    }

    // ------------------------------------------------------------------ 5
    std::printf("Test 5: flexible stem in steady current, implicit vs explicit\n");
    {
        double tip[2], Fx[2];
        for(int integ=0; integ<2; ++integ)
        {
            rodtree rt;
            load(rt,cantilever(L,r,E,1030.0,20,V3(0,0,1),0.0));
            rt.set_fluid_density(1000.0);
            rt.set_gravity(V3(0,0,-9.81));
            rt.set_integrator(integ);
            set_fluid(rt,V3(0.3,0,0),V3::Zero(),1.0);
            for(int s=0; s<4000; ++s) {rt.advance(2.0e-3); rt.compute_hydro();}
            tip[integ] = rt.tip_displacement(0)(0);
            Fx[integ] = rt.base_force(0)(0);
        }
        const double F_rigid = 0.5*1000.0*1.2*2.0*r*L*0.09;
        std::printf("    tip dx: implicit %.5f m  explicit %.5f m\n", tip[0], tip[1]);
        std::printf("    base Fx: implicit %.5f N  explicit %.5f N  rigid-stem drag %.5f N  (reconfiguration factor %.3f)\n",
                    Fx[0], Fx[1], F_rigid, Fx[0]/F_rigid);
        check(std::fabs(tip[0]-tip[1])<1.0e-3*L,"steady tip deflection agrees (1e-3 L)");
        check(std::fabs(Fx[0]-Fx[1])<1.0e-3*F_rigid,"steady base force agrees");
        check(Fx[0]<F_rigid,"flexible stem carries less drag than rigid (reconfiguration)");
    }

    // ------------------------------------------------------------------ 6
    std::printf("Test 6: branching fan in oscillatory flow, implicit vs explicit\n");
    {
        const double T = 2.0, U0 = 0.2, dt = 1.0e-3;
        std::vector<double> xs[2];
        for(int integ=0; integ<2; ++integ)
        {
            rodtree rt;
            load(rt,colony_fan(0.0,4));
            rt.set_fluid_density(1000.0);
            rt.set_integrator(integ);
            rt.set_substeps(2);
            double t = 0.0;
            for(int s=0; s<4000; ++s)
            {
                double u = U0*std::sin(2.0*PI*t/T), a = U0*2.0*PI/T*std::cos(2.0*PI*t/T);
                set_fluid(rt,V3(u,0,0),V3(a,0,0),1.0);
                rt.advance(dt); rt.compute_hydro(); t += dt;
                if(s%20==0) xs[integ].push_back(rt.tip_displacement(0)(0));
            }
        }
        double emax = 0.0, amax = 0.0;
        for(size_t i=0; i<xs[0].size(); ++i) {emax = std::max(emax,std::fabs(xs[0][i]-xs[1][i])); amax = std::max(amax,std::fabs(xs[1][i]));}
        std::printf("    max tip excursion %.4f m, max implicit-explicit difference %.2e m (%.2f%%)\n", amax, emax, 100.0*emax/amax);
        check(emax<0.02*amax,"implicit and explicit agree within 2% of amplitude");
    }

    // ------------------------------------------------------------------ 7
    std::printf("Test 7: large deflection, horizontal cantilever under heavy self-weight\n");
    {
        // uniformly loaded elastica; reference from a converged shooting
        // solution of EI theta'' + q (L - s) cos(theta) = 0 (computed below)
        const double gz = -2.7;
        const double q = rho*A*std::fabs(gz);
        // shooting on theta(0)=0, theta'(L)=0
        auto shoot = [&](double k0, double& tipz)
        {
            const int n = 20000; double h = L/n, th = 0.0, k = k0, z = 0.0;
            for(int i=0; i<n; ++i)
            {
                double s = i*h;
                // RK2 on (theta, kappa, z)
                double dth1 = k, dk1 = q*(L-s)*std::cos(th)/EI, dz1 = std::sin(th);
                double thm = th+0.5*h*dth1, km = k+0.5*h*dk1;
                double dth2 = km, dk2 = q*(L-s-0.5*h)*std::cos(thm)/EI, dz2 = std::sin(thm);
                th += h*dth2; k += h*dk2; z += h*dz2;
            }
            tipz = z; return k;
        };
        double lo = -100.0, hi = 0.0, tz = 0.0;
        for(int it=0; it<200; ++it) {double m = 0.5*(lo+hi); if(shoot(m,tz)>0.0) hi = m; else lo = m;}
        shoot(0.5*(lo+hi),tz);
        rodtree rt;
        load(rt,cantilever(L,r,E,rho,40,V3(1,0,0),0.4));
        rt.set_gravity(V3(0,0,gz));
        rt.set_fluid_density(0.0);
        for(int s=0; s<8000; ++s) rt.advance(1.0e-3);
        double w = rt.tip_displacement(0)(2);
        std::printf("    tip dz %.5f m  elastica %.5f m  (linear theory %.5f m)  rel.err %.3e\n", w, tz, -q*L*L*L*L/(8.0*EI), std::fabs(w-tz)/std::fabs(tz));
        check(std::fabs(w-tz)<0.01*std::fabs(tz),"large deflection within 1% of elastica");
    }

    // ------------------------------------------------------------------ 8
    std::printf("Test 8: very soft colony in fast current, large time step (adaptive sub-stepping)\n");
    {
        const char* bush =
            "colony bush\nmaterial 1.0e6 0.4 1050\nhydro 1.2 0.01 1.0\ndamping 0.02\nrefine 4\n"
            "node 1 0.60 0.30 0.000 0.004\nnode 2 0.60 0.30 0.080 0.0035\nnode 3 0.64 0.30 0.140 0.0025\n"
            "node 4 0.56 0.30 0.140 0.0025\nnode 5 0.60 0.34 0.150 0.0025\n"
            "edge 1 2\nedge 2 3\nedge 2 4\nedge 2 5\nclamp 1\nend\n";
        double tipx[2];
        const double dts[2] = {0.03, 0.001};
        int nmax = 0;
        for(int run=0; run<2; ++run)
        {
            rodtree rt;
            load(rt,bush);
            rt.set_fluid_density(998.2);
            const int ns = (int)std::lround(1.5/dts[run]);
            for(int s=0; s<ns; ++s)
            {
                set_fluid(rt,V3(0.5,0,0),V3::Zero(),1.0);
                rt.compute_hydro();
                rt.advance(dts[run]);
                if(run==0) nmax = std::max(nmax,rt.substeps_used());
            }
            tipx[run] = rt.tip_displacement(0)(0);
        }
        std::printf("    tip dx: dt=0.03 %.5f m  dt=0.001 %.5f m  (max sub-steps per 0.03 s step: %d)\n", tipx[0], tipx[1], nmax);
        check(std::isfinite(tipx[0]) && std::fabs(tipx[0]-tipx[1])<1.0e-3,"large-step steady state matches small-step solution");
    }

    // ------------------------------------------------------------------ 9
    std::printf("Test 9: polyp layer drag and flow-induced retraction (rigid stem)\n");
    {
        auto stem = [&](const char* extra)
        {
            std::ostringstream s;
            s<<"colony p\nmaterial 1.0e10 0.4 1000\nhydro 1.2 0.0 1.0\ndamping 0\nrefine 10\n"<<extra
             <<"node 1 0 0 0 0.004\nnode 2 0 0 0.3 0.004\nedge 1 2\nclamp 1\nend\n";
            return s.str();
        };
        auto run = [&](const std::string& in, double U, double tend, double& ext_out, std::vector<double>* ext_t=nullptr)
        {
            rodtree rt;
            load(rt,in);
            rt.set_fluid_density(1000.0);
            rt.set_gravity(V3::Zero());
            const double dt = 0.01;
            for(int s=0; s<(int)std::lround(tend/dt); ++s)
            {
                set_fluid(rt,V3(U,0,0),V3::Zero(),1.0);
                rt.compute_hydro();
                rt.advance(dt);
                if(ext_t) ext_t->push_back(rt.polyp_extension(0));
            }
            set_fluid(rt,V3(U,0,0),V3::Zero(),1.0);
            rt.compute_hydro();
            ext_out = rt.polyp_extension(0);
            return rt.hydro_force(0)(0);
        };
        double ex;
        const double F0 = run(stem(""),0.2,2.0,ex);
        const double F1 = run(stem("polyps 0.003 0.5 1.0\n"),0.2,2.0,ex);
        const double ratio_ref = (1.2*0.008 + 1.0*2.0*0.003*0.5)/(1.2*0.008);
        std::printf("    extended polyps: drag ratio %.5f  expected %.5f\n", F1/F0, ratio_ref);
        check(std::fabs(F1/F0-ratio_ref)<1.0e-4,"extended polyp layer adds Cd_p*2*h_p*phi_p of frontal width");

        std::vector<double> et;
        const double F2 = run(stem("polyps 0.003 0.5 1.0\npolyp_response 0.1 0.05 2.0\n"),0.2,20.0,ex,&et);
        const double e1_ref = std::pow(1.0/(1.0+0.01/2.0),100.0);
        std::printf("    retraction above U_r+dU (tau 2 s): ext(1 s) %.4f  expected %.4f (exp(-t/tau) %.4f);  ext(20 s) %.2e, drag/no-polyp %.5f\n",
                    et[99], e1_ref, std::exp(-0.5), ex, F2/F0);
        check(std::fabs(et[99]-e1_ref)<1.0e-6 && std::fabs(F2/F0-1.0)<1.0e-3,"polyps retract with time constant tau, drag returns to bare stem");

        const double F3 = run(stem("polyps 0.003 0.5 1.0\npolyp_response 0.1 0.05 0.0\n"),0.125,1.0,ex);
        std::printf("    half-way in the retraction ramp: ext %.3f, drag ratio %.5f  expected %.5f\n", ex, F3/(F0*0.125*0.125/0.04), 1.0+0.5*(ratio_ref-1.0));
        check(std::fabs(ex-0.5)<1.0e-5 && std::fabs(F3/(F0*0.125*0.125/0.04)-(1.0+0.5*(ratio_ref-1.0)))<1.0e-4,"partial extension scales the polyp drag");
    }

    std::printf("\n%s (%d failure%s)\n", failures==0 ? "ALL TESTS PASSED" : "TESTS FAILED", failures, failures==1?"":"s");
    return failures==0 ? 0 : 1;
}
