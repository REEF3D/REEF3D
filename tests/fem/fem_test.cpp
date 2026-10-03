// Standalone verification of the REEF3D FEM solid solver (no MPI, no REEF3D).
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         fem_test.cpp ../../src/fem_solid*.cpp -o fem_test
// Run:    ./fem_test [test]      tests: cantilever freq rotation j2 crackband drop collapse
#include"fem_solid.h"
#include<iostream>
#include<sstream>
#include<cmath>
#include<string>
#include<vector>
#include<cstdio>

typedef fem_solid::Vec3 Vec3;
static int nfail = 0;
static void check(bool ok,const std::string& what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

static fem_solid::material elastic(int id,double rho,double E,double nu)
{
    fem_solid::material m; m.id=id; m.type=fem_solid::MAT_ELASTIC; m.rho=rho; m.E=E; m.nu=nu; return m;
}

// cantilever under self weight: static tip deflection vs Timoshenko beam
static double cantilever(int full,double hg,double h,double* period=nullptr)
{
    const double L=1.0, b=0.1, t=0.1, rho=1000.0, E=1.0e8, nu=0.0, g=9.81;
    fem_solid s;
    s.set_lattice(0,0,0,h,h,h);
    s.add_material(elastic(1,rho,E,nu));
    s.add_box(0,L,0,b,0,t,1);
    s.add_fix(-1e-6,1e-6,-1,1,-1,1,true,true,true);
    s.set_element_type(full);
    s.set_hourglass(hg);
    s.set_gravity(Vec3(0,0,-g));
    s.build();

    int tip=-1;
    for(int i=0;i<s.nnode();++i)
    if(std::fabs(s.ref_pos(i)(0)-L)<1e-9 && std::fabs(s.ref_pos(i)(2)-t/2)<0.51*h && std::fabs(s.ref_pos(i)(1)-b/2)<0.51*h) {tip=i;break;}
    if(tip<0) for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-L)<1e-9) {tip=i;break;}

    const double I=b*t*t*t/12.0, q=rho*g*b*t, A=b*t, G=E/(2*(1+nu));
    const double wEB = q*L*L*L*L/(8*E*I);
    const double wT = wEB + q*L*L/(2*(5.0/6.0)*G*A);

    double dt=1.0e-3;
    if(period)
    {
        // undamped step response: the oscillation period is the natural period
        std::vector<double> tt,ww;
        for(int n=0;n<600;++n) {s.advance(dt); tt.push_back(s.time()); ww.push_back(s.pos(tip)(2)-s.ref_pos(tip)(2));}
        // mean = static deflection, period from upward crossings of the mean
        double mean=0; for(double w:ww) mean+=w; mean/=ww.size();
        std::vector<double> cr;
        for(size_t k=1;k<ww.size();++k) if(ww[k-1]<mean && ww[k]>=mean) cr.push_back(tt[k-1]+(mean-ww[k-1])/(ww[k]-ww[k-1])*dt);
        *period = cr.size()>=2 ? (cr.back()-cr.front())/(cr.size()-1) : -1;
        return -mean/wT;
    }

    s.set_damping(64.0);   // ~critical for the first mode
    for(int n=0;n<500;++n) s.advance(dt);
    const double w = -(s.pos(tip)(2)-s.ref_pos(tip)(2));
    std::printf("    %s hg=%.2f h=%.4f: tip %.5e  Timoshenko %.5e  ratio %.3f  (substeps/step %d, KE %.2e)\n",
        full?"full   ":"reduced",hg,h,w,wT,w/wT,s.last_substeps(),s.kinetic_energy());
    // support force must equal the weight
    Vec3 R=s.support_force();
    std::printf("    support force z %.4f N, weight %.4f N\n",R(2),-rho*g*L*b*t);
    return w/wT;
}

static void test_cantilever()
{
    std::cout<<"cantilever under self weight (L/t = 10)"<<std::endl;
    double r1=cantilever(1,0.0,0.025);
    double r2=cantilever(0,0.1,0.025);
    double r3=cantilever(0,0.05,0.025);
    double r4=cantilever(0,0.1,0.0125);
    check(r1>0.85 && r1<1.05,"full integration within 15% (4 elements over thickness)");
    check(r2>0.9 && r2<1.1,"reduced+hourglass within 10%");
    check(r4>0.95 && r4<1.05,"reduced+hourglass, refined, within 5%");
    (void)r3;
}

static void test_freq()
{
    std::cout<<"cantilever natural frequency"<<std::endl;
    const double L=1.0,b=0.1,t=0.1,rho=1000,E=1e8;
    const double f1 = 1.875104*1.875104/(2*M_PI)*std::sqrt(E*b*t*t*t/12/(rho*b*t*L*L*L*L));
    for(int full=0; full<2; ++full)
    {
        double T;
        double r=cantilever(full,full?0.0:0.1,0.025,&T);
        std::printf("    %s: f = %.3f Hz, Euler-Bernoulli %.3f Hz, ratio %.3f (mean defl. ratio %.3f)\n",full?"full":"reduced",1/T,f1,1/T/f1,r);
        check(std::fabs(1/T/f1-1)<0.1,std::string(full?"full":"reduced")+": frequency within 10%");
    }
}

static void test_rotation()
{
    std::cout<<"free spinning block: objectivity and energy conservation"<<std::endl;
    for(int full=0; full<2; ++full)
    {
        fem_solid s;
        s.set_lattice(0,0,0,0.05,0.05,0.05);
        s.add_material(elastic(1,1000,1e8,0.3));
        s.add_box(-0.1,0.1,-0.1,0.1,-0.1,0.1,1);
        s.set_element_type(full);
        s.set_gravity(Vec3(0,0,0));
        s.build();
        const Vec3 w(0.3,0.2,10.0);
        for(int i=0;i<s.nnode();++i) s.set_vel(i,w.cross(s.ref_pos(i)));
        const double E0=s.kinetic_energy();
        double emax=0;
        for(int n=0;n<1000;++n) {s.advance(1e-3); emax=std::max(emax,s.strain_energy());}
        const double E1=s.kinetic_energy()+s.strain_energy();
        // centrifugal strain energy is physical and small: rho w^2 r^2/E ~ 1e-5
        std::printf("    %s: KE0 %.5f  E(1s) %.5f  rel.drift %.2e  max strain energy/KE0 %.2e\n",full?"full":"reduced",E0,E1,(E1-E0)/E0,emax/E0);
        check(std::fabs(E1-E0)/E0<1e-2,"energy conserved within 1% over 1.6 revolutions");
        check(emax/E0<1e-3,"no spurious strain from large rotation");
    }
}

static void test_j2()
{
    std::cout<<"J2 plasticity, quasi-static uniaxial tension"<<std::endl;
    fem_solid s;
    const double h=0.05,L=1.0,A=0.01;
    s.set_lattice(0,0,0,h,h,h);
    fem_solid::material m; m.id=1; m.type=fem_solid::MAT_J2; m.rho=7850; m.E=2.1e11; m.nu=0.3; m.sigy=2.5e8; m.H=2.1e9; m.epsfail=0;
    s.add_material(m);
    s.add_box(0,L,0,0.1,0,0.1,1);
    s.add_fix(-1e-6,1e-6,-1,1,-1,1,true,false,false);
    s.add_fix(-1e-6,1e-6,-1e-6,1e-6,-1e-6,1e-6,true,true,true);
    s.add_fix(L-1e-6,L+1e-6,-1,1,-1,1,true,false,false);
    s.set_gravity(Vec3(0,0,0));
    s.set_damping(200.0);
    s.build();
    std::vector<int> endn;
    for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-L)<1e-9) endn.push_back(i);
    const double dt=2e-5, rate=0.02;   // strain rate 0.02/s ... loaded slowly vs wave transit 2e-4 s
    double sig_at_2pct=0;
    for(int n=0;n<10000;++n)
    {
        for(int i:endn){Vec3 p=s.pos(i); p(0)+=rate*L*dt; s.set_pos(i,p);}
        s.advance(dt);
    }
    const double eps = rate*1e4*dt;
    const double sig = s.support_force(-1,1e-6,-1,1,-1,1)(0)/A;
    // linear hardening: sigma = sigy + H_eff * eps_p, with eps_p = eps - sigma/E, H_eff = E H/(E+H) in terms of total plastic strain
    const double ep = eps - sig/m.E;
    const double sref = m.sigy + m.H*ep;
    sig_at_2pct=sig;
    std::printf("    strain %.4f: stress %.4e Pa, expected %.4e Pa (ratio %.4f)\n",eps,sig_at_2pct,sref,sig/sref);
    check(std::fabs(sig/sref-1)<0.03,"stress on the hardening branch within 3%");
}

static double crackband(double h,double* peak)
{
    fem_solid s;
    const double L=0.4, b=0.1;
    s.set_lattice(0,0,0,h,h,h);
    fem_solid::material c; c.id=1; c.type=fem_solid::MAT_CONCRETE; c.rho=2400; c.E=3e10; c.nu=0.0; c.ft=3e6; c.Gf=100; c.fc=3e7; c.Gc=1e4; c.derode=2.0;
    fem_solid::material w=c; w.id=2; w.ft=2.9e6;
    s.add_material(c); s.add_material(w);
    s.add_box(0,L,0,b,0,b,1);
    s.add_box(L/2-h/2-1e-9,L/2+h/2+1e-9,0,b,0,b,2);   // weak slice, one element thick
    s.add_fix(-1e-6,1e-6,-1,1,-1,1,true,true,true);
    s.add_fix(L-1e-6,L+1e-6,-1,1,-1,1,true,true,true);
    s.set_gravity(Vec3(0,0,0));
    s.set_damping(2000.0);
    s.build();
    std::vector<int> endn;
    for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-L)<1e-9) endn.push_back(i);
    const double dt=2e-6, du=4e-9;   // displacement per step
    double pk=0;
    for(int n=0;n<60000;++n)
    {
        for(int i:endn){Vec3 p=s.pos(i); p(0)+=du; s.set_pos(i,p);}
        s.advance(dt);
        pk=std::max(pk,s.support_force(-1,1e-6,-1,1,-1,1)(0)/(b*b));
    }
    *peak=pk;
    const double res = s.support_force(-1,1e-6,-1,1,-1,1)(0)/(b*b);
    std::printf("    h=%.4f: dissipated %.3f J, Gf*A = %.3f J, peak stress %.3e, residual %.2e\n",h,s.dissipated_energy(),100*b*b,pk,res);
    return s.dissipated_energy()/(100*b*b);
}

static void test_crackband()
{
    std::cout<<"concrete crack band: dissipated energy independent of the element size"<<std::endl;
    double p1,p2;
    double r1=crackband(0.025,&p1);
    double r2=crackband(0.0125,&p2);
    check(r1>0.85 && r1<1.2 && r2>0.85 && r2<1.2,"dissipation = Gf*A within 20% for both meshes");
    check(std::fabs(r1-r2)<0.15,"mesh objective (difference < 15%)");
    check(std::fabs(p1/2.9e6-1)<0.1,"peak stress = ft of the weak band");
}

static void test_drop()
{
    std::cout<<"block dropped onto the ground"<<std::endl;
    fem_solid s;
    s.set_lattice(0,0,0,0.05,0.05,0.05);
    s.add_material(elastic(1,2400,3e9,0.2));
    s.add_box(0,0.2,0,0.2,0.5,0.7,1);
    s.set_ground(0.0,1.0,0.5);
    s.set_damping(2.0);   // light material damping, the elastic block rings otherwise
    s.build();
    double zmin=1, emax=0;
    const double m=0.2*0.2*0.2*2400, E0=m*9.81*0.6;
    for(int n=0;n<3000;++n)
    {
        s.advance(1e-3);
        double zc=0; for(int i=0;i<s.nnode();++i) {zmin=std::min(zmin,s.pos(i)(2)); zc+=s.pos(i)(2);}
        zc/=s.nnode();
        emax=std::max(emax,s.kinetic_energy()+s.strain_energy()+m*9.81*zc);
    }
    double zc=0; for(int i=0;i<s.nnode();++i) zc+=s.pos(i)(2); zc/=s.nnode();
    std::printf("    after 3 s: mean z %.4f (resting 0.1), min z %.2e, KE %.3e, max total energy/initial %.4f\n",zc,zmin,s.kinetic_energy(),emax/E0);
    check(emax/E0<1.01,"no energy gain in contact");
    check(std::fabs(zc-0.1)<0.01,"block rests on the ground");
    check(zmin>-0.01,"penetration below 1 cm");
}

static void test_collapse()
{
    std::cout<<"weak concrete cantilever collapsing under self weight"<<std::endl;
    std::istringstream in(
        "lattice 0.05 0.05 0.05\n"
        "material 1 concrete 2400 3e9 0.2 2e4 50 2e6 5e3\n"
        "material 2 elastic 2400 3e9 0.2\n"
        "box 0 0.2 0 0.2 0 1.5 2          # column\n"
        "box 0.2 2.0 0 0.2 1.3 1.5 1      # cantilever arm\n"
        "fix -1 1 -1 1 -0.001 0.001 xyz\n"
        "ground 0.0 1.0 0.5\n"
        "contact on 1.0 0.5\n"
        "element full\n"
        "monitor tip 2.0 0.1 1.4\n");
    fem_solid s;
    s.read(in);
    s.build();
    s.info(std::cout);
    int n0=s.n_alive();
    double tmax=3.0, dt=2e-3;
    bool finite=true;
    for(int n=0;n<(int)(tmax/dt);++n)
    {
        s.advance(dt);
        if(n%250==0)
        {
            std::printf("    t %.2f  alive %d  debris %zu  KE %.3e  maxvm %.2e  tip z %.3f\n",s.time(),s.n_alive(),s.debris().size(),s.kinetic_energy(),s.max_vonmises(),s.pos(s.monitors()[0].node)(2));
        }
        for(int i=0;i<s.nnode();++i) if(!std::isfinite(s.pos(i)(2))) finite=false;
        if(!finite) break;
    }
    s.write_vtu("collapse_final.vtu");
    double zmin=1e9; for(int i=0;i<s.nnode();++i) zmin=std::min(zmin,s.pos(i)(2));
    std::printf("    final: eroded %d of %d, debris %zu, min z %.3f, KE %.3e\n",n0-s.n_alive(),n0,s.debris().size(),zmin,s.kinetic_energy());
    check(finite,"no NaN");
    check(s.n_alive()<n0,"the arm fails");
    check(zmin>-0.05,"fragments stay above the ground");
}

int main(int argc,char** argv)
{
    std::string w = argc>1 ? argv[1] : "all";
    if(w=="all"||w=="cantilever") test_cantilever();
    if(w=="all"||w=="freq") test_freq();
    if(w=="all"||w=="rotation") test_rotation();
    if(w=="all"||w=="j2") test_j2();
    if(w=="all"||w=="crackband") test_crackband();
    if(w=="all"||w=="drop") test_drop();
    if(w=="all"||w=="collapse") test_collapse();
    std::cout<<(nfail ? "FAILED: " : "all tests passed")<<(nfail? std::to_string(nfail):"")<<std::endl;
    return nfail ? 1 : 0;
}
