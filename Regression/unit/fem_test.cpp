// Standalone verification of the REEF3D FEM solid solver (no MPI, no REEF3D).
// Architect: Hans Bihs
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         fem_test.cpp ../../src/fem_solid*.cpp -o fem_test          (add -fopenmp for the threads test)
// Run:    ./fem_test [test]      tests: cantilever freq rotation j2 crackband drop collapse snap patch
//                                presets settle snapbeam damping rigid walls impact rebar tie rc threads (default: all)
#include"fem_solid.h"
#include<iostream>
#include<sstream>
#include<cmath>
#include<string>
#include<vector>
#include<cstdio>
#include<cstring>

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
    std::printf("    final: eroded %d of %d, debris %zu, min z %.3f, KE %.3e, rigid fragments %d, substeps per step %d\n",n0-s.n_alive(),n0,s.debris().size(),zmin,s.kinetic_energy(),s.n_rigid(),s.last_substeps());
    check(finite,"no NaN");
    check(s.n_alive()<n0,"the arm fails");
    check(zmin>-0.05,"fragments stay above the ground");
    check(s.n_rigid()>=1 && s.rigid(0).fragment,"the broken arm becomes a rigid body (fragments rigid)");
    check(s.kinetic_energy()<1.0 && zmin>-0.005,"the fragment comes to rest on the ground (penetration below 5 mm)");
}


// ----------------------------------------------------------------------
// accessibility layer: snapping, presets, settle, check
// ----------------------------------------------------------------------
#include<fstream>
static void write_cylinder_stl(const std::string& fn,double xc,double yc,double r,double z0,double z1,int n)
{
    std::ofstream f(fn.c_str());
    f<<"solid c\n";
    auto tri=[&](double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz){
        f<<"facet normal 0 0 0\nouter loop\nvertex "<<ax<<" "<<ay<<" "<<az<<"\nvertex "<<bx<<" "<<by<<" "<<bz<<"\nvertex "<<cx<<" "<<cy<<" "<<cz<<"\nendloop\nendfacet\n";};
    for(int i=0;i<n;++i)
    {
        double a0=2*M_PI*i/n, a1=2*M_PI*(i+1)/n;
        double x0=xc+r*cos(a0), y0=yc+r*sin(a0), x1=xc+r*cos(a1), y1=yc+r*sin(a1);
        tri(x0,y0,z0,x1,y1,z0,x1,y1,z1); tri(x0,y0,z0,x1,y1,z1,x0,y0,z1);
        tri(xc,yc,z0,x1,y1,z0,x0,y0,z0); tri(xc,yc,z1,x0,y0,z1,x1,y1,z1);
    }
    f<<"endsolid c\n";
}

static void test_snap()
{
    std::cout<<"surface snapping: cylinder r = 0.1 m, h = 0.2 m, lattice 0.02 m"<<std::endl;
    write_cylinder_stl("cyl.stl",0,0,0.1,0,0.2,256);
    const double Vex = M_PI*0.01*0.2*(256.0/(2*M_PI)*sin(2*M_PI/256.0));   // polygon volume
    double err[2], aerr[2];
    for(int snap=0; snap<2; ++snap)
    {
        std::istringstream in("lattice 0.02 0.02 0.02 0.001 0.001 0\nmaterial concrete C30\nstl cyl.stl\nfix base\nsnap "+std::string(snap?"on":"off")+"\n");
        fem_solid s; s.read(in); s.build();
        fem_solid::check_info ci=s.check();
        err[snap]=ci.volume/Vex-1;
        aerr[snap]=ci.surface_area/(2*M_PI*0.1*0.2+2*M_PI*0.01)-1;
        std::printf("    snap %-3s: volume %.6f m3 (exact %.6f), error %+.2f %%, surface %.4f m2 (exact %.4f)\n",snap?"on":"off",ci.volume,Vex,100*err[snap],ci.surface_area,2*M_PI*0.1*0.2+2*M_PI*0.01);
    }
    check(std::fabs(err[1])<0.01,"snapped volume within 1 %");
    check(std::fabs(aerr[1])<0.02 && std::fabs(aerr[0])>0.1,"snapping removes the stair steps (surface area within 2 %, voxel +18 %)");
}

static void test_patch()
{
    std::cout<<"patch test on the snapped mesh: uniform strain gives no interior forces"<<std::endl;
    std::istringstream in("lattice 0.02 0.02 0.02 0.001 0.001 0\nmaterial 1 elastic 2000 1e9 0.25\nstl cyl.stl\n");
    fem_solid s; s.read(in); s.set_gravity(Vec3(0,0,0)); s.set_contact(false,1,0); s.build();
    // linear displacement field
    for(int i=0;i<s.nnode();++i){Vec3 X=s.ref_pos(i); s.set_pos(i,X+Vec3(1e-4*X(0)+2e-4*X(1),-1e-4*X(2),3e-4*X(0)));}
    const double dt=1e-7; s.advance(dt);
    std::vector<char> surf(s.nnode(),0);
    for(const auto& f: s.surface()) for(int q=0;q<4;++q) surf[f.n[q]]=1;
    double vin=0,vsurf=0;
    for(int i=0;i<s.nnode();++i){double vv=s.vel(i).norm(); if(surf[i]) vsurf=std::max(vsurf,vv); else vin=std::max(vin,vv);}
    std::printf("    max |v| interior %.3e, surface %.3e (ratio %.1e)\n",vin,vsurf,vin/vsurf);
    check(vin<1e-6*vsurf,"interior nodes in equilibrium (patch test passed)");
}

static void test_presets()
{
    std::cout<<"material presets"<<std::endl;
    fem_solid::material m;
    bool ok=fem_solid::preset("concrete","C30/37",m);
    std::printf("    concrete C30/37: E %.3g ft %.3g fc %.3g Gf %.1f Gc %.0f rho %.0f\n",m.E,m.ft,m.fc,m.Gf,m.Gc,m.rho);
    check(ok && m.E==33e9 && m.ft==2.9e6 && m.fc==38e6 && std::fabs(m.Gf-140.8)<1.0,"C30/37 values (EN 1992, MC2010)");
    ok=fem_solid::preset("steel","S355",m);
    check(ok && m.sigy==355e6 && m.type==fem_solid::MAT_J2,"S355");
    check(!fem_solid::preset("concrete","C99",m),"unknown preset rejected");
    std::istringstream in("material concrete C25 ft 2.0e6\nmaterial steel S235\nlattice 0.1 0.1 0.1\nbox 0 1 0 1 0 1 1\nbox 1 2 0 1 0 1\n");
    fem_solid s; s.read(in); s.build();
    check(s.material_count()==2 && s.mat(0).ft==2.0e6 && s.mat(1).id==2,"auto ids, override, default material for shapes");
}

static void test_settle_check()
{
    std::cout<<"settle and check: cantilever under self weight (L/t = 10)"<<std::endl;
    std::istringstream in("lattice 0.025 0.025 0.025\nmaterial 1 elastic 1000 1e8 0\nbox 0 1 0 0.1 0 0.1\nfix -1e-6 1e-6 -1 1 -1 1\nmonitor tip auto\n");
    fem_solid s; s.read(in); s.build();
    double res; bool ok=s.settle(200000,1e-6,&res);
    const double L=1,b=0.1,t=0.1,rho=1000,E=1e8,g=9.81, I=b*t*t*t/12, q=rho*g*b*t, G=E/2;
    const double wT=q*L*L*L*L/(8*E*I)+q*L*L/(2*(5.0/6.0)*G*b*t);
    int tip=-1; for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-1)<1e-9 && std::fabs(s.ref_pos(i)(2)-0.05)<1e-9 && std::fabs(s.ref_pos(i)(1)-0.05)<1e-9) tip=i;
    const double w=-(s.pos(tip)(2)-s.ref_pos(tip)(2));
    std::printf("    settle %s (residual %.1e): tip %.5e m, Timoshenko %.5e, ratio %.4f\n",ok?"converged":"NOT converged",res,w,wT,w/wT);
    check(ok && std::fabs(w/wT-1)<0.02,"settled deflection within 2 %");
    std::istringstream in2("lattice 0.025 0.025 0.025\nmaterial 1 elastic 1000 1e8 0\nbox 0 1 0 0.1 0 0.1\nfix -1e-6 1e-6 -1 1 -1 1\n");
    fem_solid s2; s2.read(in2); s2.build();
    fem_solid::check_info ci=s2.check();
    s2.write_check(std::cout,ci);
    const double f1 = 1.875104*1.875104/(2*M_PI)*std::sqrt(E*I/(rho*b*t*L*L*L*L));
    std::printf("    Rayleigh f_z %.3f Hz, Euler-Bernoulli %.3f Hz, ratio %.3f\n",ci.freq[2],f1,ci.freq[2]/f1);
    check(ci.freq_ok[2] && std::fabs(ci.freq[2]/f1-1)<0.03,"first frequency within 3 %");
    check(std::fabs(ci.sw_support(2)+rho*g*L*b*t)<1e-3*rho*g*L*b*t,"support force = weight");
}

static void test_snapped_cantilever()
{
    std::cout<<"cantilever with a thickness off the lattice (t = 0.09 m, h = 0.025 m, snapped)"<<std::endl;
    std::istringstream in("lattice 0.025 0.025 0.025\nmaterial 1 elastic 1000 1e8 0\nbox 0 1 0 0.1 0 0.09\nfix -1e-6 1e-6 -1 1 -1 1\n");
    fem_solid s; s.read(in); s.build();
    double res; s.settle(400000,1e-6,&res);
    const double L=1,b=0.1,t=0.09,rho=1000,E=1e8,g=9.81, I=b*t*t*t/12, q=rho*g*b*t, G=E/2;
    const double wT=q*L*L*L*L/(8*E*I)+q*L*L/(2*(5.0/6.0)*G*b*t);
    double w=0; for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-1)<1e-9) w=std::max(w,-(s.pos(i)(2)-s.ref_pos(i)(2)));
    fem_solid::check_info ci=s.check();
    std::printf("    volume %.5f (exact %.5f), tip %.5e, Timoshenko %.5e, ratio %.3f\n",ci.volume,L*b*t,w,wT,w/wT);
    check(std::fabs(ci.volume/(L*b*t)-1)<1e-6,"exact volume after snapping");
    check(std::fabs(w/wT-1)<0.06,"deflection within 6 % with distorted elements");
}

static void test_damping()
{
    std::cout<<"structural damping: 5 % at the first frequency, free parts undamped"<<std::endl;
    const double L=1.0,b=0.1,t=0.1,rho=1000,E=1e8,g=9.81,h=0.025;
    fem_solid s;
    s.set_lattice(0,0,0,h,h,h);
    s.add_material(elastic(1,rho,E,0.0));
    s.add_box(0,L,0,b,0,t,1);                       // cantilever
    s.add_box(2.0,2.2,0,0.2,0,0.2,1);               // free block, separate body
    s.add_fix(-1e-6,1e-6,-1,1,-1,1,true,true,true);
    s.set_gravity(Vec3(0,0,-g));
    s.set_damping_ratio(0.05);
    s.build();
    s.prepare_damping();
    const double f1 = 1.875104*1.875104/(2*M_PI)*std::sqrt(E*t*t/12.0/(rho*L*L*L*L));
    std::printf("    damping frequency %.3f Hz (Rayleigh), Euler-Bernoulli %.3f Hz, alpha %.3f 1/s\n",s.damping_frequency(),f1,s.damping_alpha());
    check(std::fabs(s.damping_frequency()/f1-1)<0.05,"Rayleigh frequency of the supported part within 5%");

    int tip=-1;
    for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-L)<1e-9) {tip=i;break;}
    std::vector<int> blk;
    for(int i=0;i<s.nnode();++i) if(s.ref_pos(i)(0)>1.5) blk.push_back(i);
    const Vec3 c(2.1,0.1,0.1), w(0,0,5.0);
    for(int i:blk) s.set_vel(i,w.cross(s.ref_pos(i)-c));
    auto ke_rot=[&](){double e=0; Vec3 vm=Vec3::Zero(); double M=0; for(int i:blk){vm+=s.mass(i)*s.vel(i); M+=s.mass(i);} vm/=M;
                      for(int i:blk) e+=0.5*s.mass(i)*(s.vel(i)-vm).squaredNorm(); return e;};
    const double E0=ke_rot();

    const double dt=1e-3;
    std::vector<double> ww;
    for(int n=0;n<1500;++n) {s.advance(dt); ww.push_back(s.pos(tip)(2)-s.ref_pos(tip)(2));}
    // peaks of the deflection about the final static value
    const double ws = ww.back();
    std::vector<double> pk;
    for(size_t k=1;k+1<ww.size();++k) if(ww[k]-ws<ww[k-1]-ws && ww[k]-ws<=ww[k+1]-ws && ws-ww[k]>0) pk.push_back(ws-ww[k]);
    double zm=-1;
    if(pk.size()>=3)
    {
        const double delta = std::log(pk[0]/pk[2])/2.0;
        zm = delta/std::sqrt(4*M_PI*M_PI+delta*delta);
    }
    Vec3 vm=Vec3::Zero(); double M=0; for(int i:blk){vm+=s.mass(i)*s.vel(i); M+=s.mass(i);} vm/=M;
    std::printf("    measured damping ratio %.4f (log decrement of the tip), free block: v_z %.4f m/s (free fall %.4f), rotational energy %.5f -> %.5f\n",
        zm,vm(2),-g*s.time(),E0,ke_rot());
    check(zm>0.04 && zm<0.06,"damping ratio of the first mode 5 +- 1 %");
    check(std::fabs(vm(2)/(-g*s.time())-1)<1e-3,"free part falls undamped");
    check(std::fabs(ke_rot()/E0-1)<0.01,"rotation of the free part undamped");
}

static void test_rigid()
{
    std::cout<<"rigid bodies: free fall, torque-free rotation, resting on the ground, timber on a cantilever"<<std::endl;
    {
        // free fall + spin about a non-principal axis: momentum and angular momentum conserved
        fem_solid s;
        s.set_lattice(0,0,0,0.05,0.05,0.05);
        fem_solid::material m=elastic(1,500,1e10,0.3); m.rigid=true;
        s.add_material(m);
        s.add_box(0,0.4,0,0.2,0,0.1,1);
        s.set_gravity(Vec3(0,0,-9.81));
        s.build();
        check(s.n_rigid()==1,"one rigid body");
        const fem_solid::rigid_body& rb=s.rigid(0);
        const Vec3 w0(1.0,2.0,3.0);
        // initial spin: set through the angular momentum
        const_cast<fem_solid::rigid_body&>(rb).w = w0;
        const_cast<fem_solid::rigid_body&>(rb).L = rb.I0*w0;
        double E0=0.5*w0.dot(rb.I0*w0);
        for(int n=0;n<1000;++n) s.advance(1e-3);
        const Eigen::Matrix3d I=rb.R*rb.I0*rb.R.transpose();
        const double E1=0.5*rb.w.dot(I*rb.w);
        const double dz=rb.c(2)-rb.c0(2), ex=-0.5*9.81*s.time()*s.time();
        std::printf("    substeps per ms %d, free fall %.5f m (exact %.5f), rotational energy drift %.2e, |R^T R - I| %.1e\n",
            s.last_substeps(),dz,ex,(E1-E0)/E0,(rb.R.transpose()*rb.R-Eigen::Matrix3d::Identity()).norm());
        check(std::fabs(dz/ex-1)<1e-3,"free fall");
        check(std::fabs(E1/E0-1)<2e-2,"rotational energy conserved within 2 % (1 s, 3 rad/s)");
        check((rb.R.transpose()*rb.R-Eigen::Matrix3d::Identity()).norm()<1e-10,"rotation stays orthonormal");
        // shape kept: distance of two nodes
        const double d0=(s.ref_pos(0)-s.ref_pos(s.nnode()-1)).norm(), d1=(s.pos(0)-s.pos(s.nnode()-1)).norm();
        check(std::fabs(d1-d0)<1e-9,"rigid: node distances kept");
    }
    {
        // dropped on the ground: comes to rest at the ground level
        fem_solid s;
        s.set_lattice(0,0,0,0.05,0.05,0.05);
        fem_solid::material m=elastic(1,2400,3e10,0.2); m.rigid=true;
        s.add_material(m);
        s.add_box(0,0.3,0,0.3,0.2,0.4,1);
        s.set_gravity(Vec3(0,0,-9.81));
        s.set_ground(0.0,1.0,0.5);
        s.build();
        for(int n=0;n<1500;++n) s.advance(1e-3);
        const fem_solid::rigid_body& rb=s.rigid(0);
        double zmin=1e9; for(int i=0;i<s.nnode();++i) zmin=std::min(zmin,s.pos(i)(2));
        std::printf("    block on the ground: lowest node %.2e m, speed %.2e m/s, substeps per ms %d\n",zmin,rb.V.norm(),s.last_substeps());
        check(zmin>-5e-3 && zmin<1e-3,"rests on the ground (penetration below 5 mm)");
        check(rb.V.norm()<1e-2,"at rest");
    }
    {
        // rigid block falling onto an elastic cantilever: contact between a rigid and a deformable body
        fem_solid s;
        s.set_lattice(0,0,0,0.025,0.025,0.025);
        s.add_material(elastic(1,1000,1e8,0.0));
        fem_solid::material m=elastic(2,500,1e9,0.3); m.rigid=true;
        s.add_material(m);
        s.add_box(0,1.0,0,0.1,0,0.1,1);
        s.add_box(0.8,0.9,0,0.1,0.15,0.25,2);
        s.add_fix(-1e-6,1e-6,-1,1,-1,1,true,true,true);
        s.set_gravity(Vec3(0,0,-9.81));
        s.set_damping(20.0);
        s.build();
        check(s.n_rigid()==1,"block rigid, cantilever deformable");
        for(int n=0;n<1500;++n) s.advance(1e-3);
        const fem_solid::rigid_body& rb=s.rigid(0);
        int tip=-1; for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(0)-1.0)<1e-9 && s.ref_pos(i)(2)>0.09) {tip=i;break;}
        const Vec3 R=s.support_force();
        const double W=9.81*(1000*1.0*0.1*0.1+500*0.1*0.1*0.1);
        std::printf("    block on the cantilever: block centre z %.4f (start %.4f), tip deflection %.4f m, support force %.3f N, total weight %.3f N\n",
            rb.c(2),rb.c0(2),s.pos(tip)(2)-s.ref_pos(tip)(2),R(2),W);
        check(rb.c(2)<rb.c0(2)-0.04 && rb.c(2)>0.1,"block lands on the beam and stays on it");
        check(std::fabs(-R(2)/W-1)<0.05,"support carries beam + block within 5 %");
    }
}

static void test_walls()
{
    std::cout<<"contact with walls and an inclined bed (sampled level set): stick and slide"<<std::endl;
    for(int c=0;c<2;++c)
    {
        const double deg = c==0 ? 20.0 : 35.0, mu = 0.5, g = 9.81;
        const double a = deg*M_PI/180.0;
        const Vec3 n(-std::sin(a),0.0,std::cos(a));       // bed normal, bed through the origin
        fem_solid s;
        s.set_lattice(0,0,0,0.05,0.05,0.05);
        fem_solid::material m=elastic(1,2000,1e9,0.3); m.rigid=true;
        s.add_material(m);
        s.add_box(0,0.2,0,0.2,0,0.1,1);
        s.set_gravity(Vec3(0,0,-g));
        s.set_contact(true,1.0,mu);
        s.set_ground(-100.0,1.0,mu);                      // only to set the friction, far below
        s.set_bed_contact(true);
        s.build();
        // place the block on the bed: rotate the reference frame by the slope
        fem_solid::rigid_body& rb=const_cast<fem_solid::rigid_body&>(s.rigid(0));
        rb.R = Eigen::AngleAxisd(-a,Vec3(0,1,0)).toRotationMatrix();
        rb.c = Vec3(0,0.1,0) + 0.05*n;                    // resting on the bed (half thickness)
        for(int i=0;i<s.nnode();++i) s.set_pos(i, rb.c + rb.R*(s.ref_pos(i)-rb.c0));
        double t=0;
        for(int k=0;k<600;++k)
        {
            s.clear_bed_samples();
            for(int i=0;i<s.nnode();++i) s.set_bed_sample(i, n.dot(s.pos(i)), n);
            s.advance(1e-3); t+=1e-3;
        }
        const Vec3 tdir(std::cos(a),0.0,std::sin(a));     // up the slope
        const double sdist = -(rb.c-(Vec3(0,0.1,0)+0.05*n)).dot(tdir);
        const double aexp = std::max(0.0, g*(std::sin(a)-mu*std::cos(a)));
        const double sexp = 0.5*aexp*t*t;
        double pen=0; for(int i=0;i<s.nnode();++i) pen=std::max(pen,-n.dot(s.pos(i)));
        std::printf("    slope %.0f deg, mu %.1f: slid %.4f m in %.2f s (Coulomb %.4f m), max penetration %.1e m\n",deg,mu,sdist,t,sexp,pen);
        if(c==0) check(std::fabs(sdist)<5e-3,"sticks on 20 deg (tan < mu)");
        else check(std::fabs(sdist/sexp-1)<0.05,"slides on 35 deg with a = g (sin - mu cos) within 5 %");
        check(pen<5e-3,"no penetration of the bed (below 5 mm)");
    }
    {
        // thrown against a wall plane: bounces back, does not pass
        fem_solid s;
        s.set_lattice(0,0,0,0.05,0.05,0.05);
        fem_solid::material m=elastic(1,500,1e9,0.3); m.rigid=true;
        s.add_material(m);
        s.add_box(0,0.2,0,0.2,0,0.2,1);
        s.set_gravity(Vec3(0,0,0));
        s.add_contact_plane(Vec3(-1,0,0),-0.5);           // wall at x = 0.5, inside is x < 0.5
        s.build();
        fem_solid::rigid_body& rb=const_cast<fem_solid::rigid_body&>(s.rigid(0));
        rb.V = Vec3(2.0,0,0);
        double xmax=0;
        for(int k=0;k<500;++k){ s.advance(1e-3); for(int i=0;i<s.nnode();++i) xmax=std::max(xmax,s.pos(i)(0)); }
        std::printf("    wall impact at 2 m/s: max overshoot %.2e m, rebound velocity %.3f m/s, max contact force %.1f N\n",xmax-0.5,rb.V(0),rb.fcmax);
        check(xmax-0.5<0.01,"stops at the wall (overshoot below 1 cm)");
        check(rb.V(0)<0.0 && rb.V(0)>-2.0,"rebounds with less than the impact speed");
    }
}

// one rigid block (and optionally a second one or a fixed elastic wall) thrown at a
// wall: peak contact force, duration and rebound of the debris impact
struct impact_result {double Fmax, dur, vreb, k, Fsup;};
static impact_result impact_run(double k,double zeta,double cap,double A,int partner,double u=2.0)
{
    // block 0.4 x 0.2 x 0.2 m, rho 500 (8 kg), stiffness k [N/m] (k<=0: bar from E = 1e10)
    fem_solid s;
    s.set_lattice(0,0,0,0.05,0.05,0.05);
    fem_solid::material m=elastic(1,500,1e10,0.3); m.rigid=true; m.kdebris=k>0?k:0.0; m.fcrush=cap;
    s.add_material(m);
    s.add_box(0,0.4,0,0.2,0,0.2,1);
    if(partner==1)
    {
        // second rigid block 0.2 x 0.2 x 0.2 m (4 kg), twice as stiff, coming the other way
        fem_solid::material m2=m; m2.id=2; m2.kdebris=2.0*k;
        s.add_material(m2);
        s.add_box(0.45,0.65,0,0.2,0,0.2,2);
    }
    else if(partner==2)
    {
        // fixed concrete-like wall (E 30 GPa) 0.2 m thick
        s.add_material(elastic(2,2400,3e10,0.2));
        s.add_box(0.45,0.65,-0.1,0.3,-0.1,0.3,2);
        s.add_fix(0.65-1e-6,0.65+1e-6,-1,1,-1,1,true,true,true);
    }
    else
    s.add_contact_plane(Vec3(-1,0,0),-0.45);         // wall at x = 0.45
    s.set_gravity(Vec3(0,0,0));
    s.set_debris_damping(zeta);
    s.set_damping_ratio(0.0);
    s.build();
    fem_solid::rigid_body& rb=const_cast<fem_solid::rigid_body&>(s.rigid(0));
    rb.V = Vec3(u,0,0);
    if(partner==1) const_cast<fem_solid::rigid_body&>(s.rigid(1)).V = Vec3(-u,0,0);
    if(A>0.0) s.set_rigid_added_mass(0,A);
    impact_result r{0,0,0,rb.k,0};
    const double dt=2e-5;
    double t0=-1, t1=-1;
    for(int n=0;n<5000;++n)
    {
        s.advance(dt);
        const double F=rb.fcstep;
        if(F>0 && t0<0) t0=s.time()-dt;
        if(F>0) t1=s.time();
        r.Fmax=std::max(r.Fmax,F);
        if(partner==2) r.Fsup=std::max(r.Fsup,std::fabs(s.support_force()(0)));
        if(t0>0 && F==0 && s.time()>t1+0.02) break;
    }
    r.dur=t1-t0; r.vreb=rb.V(0);
    return r;
}

static void test_impact()
{
    std::cout<<"debris impact: rigid body as a spring k, peak u sqrt(k M), duration pi sqrt(M/k)"<<std::endl;
    const double M=8.0, u=2.0, k=1.0e6;
    const double F0=u*std::sqrt(k*M), T0=M_PI*std::sqrt(M/k);
    {
        impact_result r=impact_run(k,0.0,0.0,0.0,0);
        std::printf("    wall, elastic: peak %.1f N (u sqrt(kM) %.1f), duration %.2f ms (%.2f), rebound %.3f m/s\n",r.Fmax,F0,1e3*r.dur,1e3*T0,r.vreb);
        check(std::fabs(r.Fmax/F0-1)<0.02,"peak force u sqrt(k M) within 2 %");
        check(std::fabs(r.dur/T0-1)<0.03,"impact duration pi sqrt(M/k) within 3 %");
        check(std::fabs(r.vreb/u+1)<0.01,"elastic rebound -u within 1 %");
    }
    {
        impact_result r=impact_run(k,0.5,0.0,0.0,0);
        std::printf("    wall, damping 50 %% (unloading): peak %.1f N, duration %.2f ms, restitution %.3f\n",r.Fmax,1e3*r.dur,-r.vreb/u);
        check(std::fabs(r.Fmax/F0-1)<0.02,"damping on unloading only: peak unchanged within 2 %");
        check(-r.vreb/u>0.50 && -r.vreb/u<0.60,"restitution 0.55 +- 0.05");
    }
    {
        impact_result r=impact_run(k,0.0,0.0,20.0,0);
        std::printf("    wall, added mass 20 kg on the 8 kg body: peak %.1f N, rebound %.3f m/s\n",r.Fmax,r.vreb);
        check(std::fabs(r.Fmax/F0-1)<0.02,"contact acts on the body mass alone (added mass of the fluid loads not in the impact)");
    }
    {
        impact_result r=impact_run(k,0.0,0.5*F0,0.0,0);
        std::printf("    wall, crushing force %.1f N: peak %.1f N, rebound %.3f m/s\n",0.5*F0,r.Fmax,r.vreb);
        check(r.Fmax<=0.5*F0*1.001,"force capped at the crushing force");
        check(std::fabs(r.vreb/u+0.5)<0.03,"rebound with the elastic part of the energy only (u/2 for a cap at half the peak)");
    }
    {
        impact_result r=impact_run(0.0,0.0,0.0,0.0,0,1.0);
        const double kb=1e10*(0.2*0.2)/0.4;
        std::printf("    bar model E A / L: k %.3e N/m (%.3e), peak at 1 m/s %.0f N (%.0f)\n",r.k,kb,r.Fmax,std::sqrt(kb*M));
        check(std::fabs(r.k/kb-1)<1e-9,"stiffness from the axial bar E A / L");
        check(std::fabs(r.Fmax/std::sqrt(kb*M)-1)<0.03,"peak with the bar stiffness within 3 %");
    }
    {
        impact_result r=impact_run(k,0.0,0.0,0.0,1);
        const double ks=k*2*k/(3*k), me=M*4.0/12.0, F2=2*u*std::sqrt(ks*me);
        std::printf("    two bodies head-on (8 kg at 2 m/s, 4 kg at -2 m/s, k and 2k): peak %.1f N (series %.1f)\n",r.Fmax,F2);
        check(std::fabs(r.Fmax/F2-1)<0.03,"two bodies: springs in series, reduced mass, within 3 %");
    }
    {
        // tilted block dropped on the ground with a low crushing force: the crushed corner
        // must not let the other corners sink in when it rocks back (permanent set per point)
        fem_solid s;
        s.set_lattice(0,0,0,0.05,0.05,0.05);
        fem_solid::material m=elastic(1,500,1e10,0.3); m.rigid=true; m.kdebris=1e6; m.fcrush=2000.0;
        s.add_material(m);
        s.add_box(0,0.4,0,0.2,0,0.2,1);
        s.set_gravity(Vec3(0,0,-9.81));
        s.set_ground(0.0,1.0,0.5);
        s.set_debris_damping(0.95);                         // no rebound: the corner stays in contact
        s.build();
        fem_solid::rigid_body& rb=const_cast<fem_solid::rigid_body&>(s.rigid(0));
        rb.R = Eigen::AngleAxisd(0.35,Vec3(0,1,0)).toRotationMatrix();
        rb.c = Vec3(0.2,0.1,0.45);
        for(int i=0;i<s.nnode();++i) s.set_pos(i, rb.c + rb.R*(s.ref_pos(i)-rb.c0));
        double pmax=0;
        for(int n=0;n<3000;++n)
        {
            s.advance(1e-3);
            for(int i=0;i<s.nnode();++i) pmax=std::max(pmax,-s.pos(i)(2));
        }
        double zmin=1e9; for(int i=0;i<s.nnode();++i) zmin=std::min(zmin,s.pos(i)(2));
        std::printf("    tilted block (8 kg) on the ground, crushing at 2 kN: max penetration %.4f m, at rest %.2e m, speed %.2e m/s\n",pmax,-zmin,rb.V.norm());
        check(pmax>0.005,"the corner crushes");
        check(-zmin<1e-3 && rb.V.norm()<1e-2,"comes to rest flat on the ground (no sinking of the other corners)");
    }
    {
        impact_result r=impact_run(k,0.0,0.0,0.0,2);
        std::printf("    fixed concrete wall (FEM): contact %.1f N, support force %.1f N (u sqrt(kM) %.1f)\n",r.Fmax,r.Fsup,F0);
        check(std::fabs(r.Fmax/F0-1)<0.03,"debris much softer than the structure: contact u sqrt(k M) within 3 %");
        check(r.Fsup>0.8*F0 && r.Fsup<1.3*F0,"support force of the stiff wall close to the contact force");
    }
}

// ----------------------------------------------------------------------
// reinforced concrete: smeared bars
// ----------------------------------------------------------------------

static void test_rebar_input()
{
    std::cout<<"reinforcement input: presets, skin placement, steel volume"<<std::endl;
    std::istringstream in(
        "lattice 0.05 0.05 0.05\n"
        "material 1 concrete C30\n"
        "rebar B500B z 1% x 0.25%\n"
        "box 0 0.5 0 0.5 0 1.0\n"
        "fix -1 1 -1 1 -0.001 0.001 xyz\n");
    fem_solid s;
    s.read(in);
    s.build();
    s.info(std::cout);
    double vz=0, vx=0, V=0; int nskin=0;
    for(int e=0;e<s.nelem();++e)
    {
        const double ve=0.05*0.05*0.05;
        V+=ve; vz+=s.elem(e).rs[2]*ve; vx+=s.elem(e).rs[0]*ve;
        if(s.elem(e).iz==0 && s.elem(e).rs[2]>0) ++nskin;
    }
    fem_solid::material m; fem_solid::rebar_preset("B500B",m);
    std::printf("    skin elements per layer %d of 100, local z fraction %.4f, steel volume z %.4f %% x %.4f %%, B500B f_y %.0f MPa H %.0f MPa eps_u %.3f\n",
                nskin,s.elem(0).rs[2],100*vz/V,100*vx/V,m.sfy/1e6,m.sH/1e6,m.seu);
    check(nskin==36,"the bars lie in the outer element layer (36 of 100 elements)");
    check(std::fabs(vz/V-0.01)<1e-12 && std::fabs(vx/V-0.0025)<1e-12,"the steel volume of the member is the given fraction");
    check(std::fabs(m.sfy-550e6)<1 && std::fabs(m.seu-0.05)<1e-12,"B500B: f_y = 1.1 x 500 MPa, rupture at 5 %");
}

static void test_rebar_tension()
{
    std::cout<<"reinforced tie in tension: cracking, yield plateau rho f_y, rupture at eps_u"<<std::endl;
    const double h=0.05, L=0.5;
    fem_solid s;
    s.set_lattice(0,0,0,h,h,h);
    s.set_plane_strain(true);
    fem_solid::material c; c.id=1; c.type=fem_solid::MAT_CONCRETE; c.rho=2400; c.E=33e9; c.nu=0.0;
    c.ft=2.9e6; c.Gf=140; c.fc=38e6; c.Gc=54000;
    fem_solid::rebar_preset("B500B",c); c.rs[2]=0.01; c.rskin=0.0;
    s.add_material(c);
    s.add_box(0,h,0,h,0,L,1);
    s.add_fix(-1,1,-1,1,-1e-6,1e-6,true,true,true);
    s.set_gravity(Vec3(0,0,0));
    s.set_damping(200.0);
    s.build();
    std::vector<int> top;
    for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(2)-L)<1e-9) top.push_back(i);
    const double dt=4e-6, v=0.02;
    double u=0, s2=0, s4=0, pk=0;
    while(u<0.03)
    {
        u+=v*dt;
        for(int i:top){Vec3 p=s.ref_pos(i); p(2)+=u; s.set_pos(i,p); s.set_vel(i,Vec3(0,0,v));}
        s.advance(dt);
        const double sg=s.support_force(-1,1,-1,1,-1e-6,1e-6)(2)/(h*h);
        pk=std::max(pk,sg);
        if(std::fabs(u/L-0.02)<0.5*v*dt/L) s2=sg;
        if(std::fabs(u/L-0.04)<0.5*v*dt/L) s4=sg;
    }
    const double sf=s.support_force(-1,1,-1,1,-1e-6,1e-6)(2)/(h*h);
    const double fy=c.sfy, H=c.sH, Es=c.sE;
    // linear kinematic hardening on the Green-Lagrange strain E = e + e^2/2: S = f_y + H E_p,
    // E_p = (E - f_y/E_s)/(1 + H/E_s); nominal stress (1 + e) S
    auto bar=[&](double e){const double G=e+0.5*e*e; return (1+e)*0.01*(fy + H*(G-fy/Es)/(1+H/Es));};
    std::printf("    stress at 2 %% strain %.3f MPa (bar %.3f), at 4 %% %.3f (bar %.3f), peak %.3f, after 6 %% %.2e MPa, ruptured %d\n",s2/1e6,bar(0.02)/1e6,s4/1e6,bar(0.04)/1e6,pk/1e6,sf/1e6,s.bars_ruptured());
    check(std::fabs(s2/bar(0.02)-1)<0.03 && std::fabs(s4/bar(0.04)-1)<0.03,"cracked tie: the bars carry rho sigma_s (within 3 %)");
    check(s.bars_ruptured()>0 && std::fabs(sf)<1e-3*pk,"the bars rupture beyond eps_u, the tie carries nothing");
}

static void rc_push(double h,int full,double ratio,double umax,double v,double* Mcr,double* My,double* uy,double* M3,double* Mpk,double* Mend,int* rupt,int* eroded)
{
    // plane-strain cantilever 0.5 x 3 m pushed at the top through a stiff spring
    const double L=3.0, D=0.5, B=h;
    fem_solid s;
    s.set_lattice(0,0,0,h,h,h);
    s.set_plane_strain(true);
    fem_solid::material c; c.id=1; c.type=fem_solid::MAT_CONCRETE; c.rho=2400; c.E=33e9; c.nu=0.0;
    c.ft=2.9e6; c.Gf=73*std::pow(38.0,0.18); c.fc=38e6; c.Gc=8.8*std::sqrt(38.0)*1000;
    if(ratio>0) {fem_solid::rebar_preset("B500B",c); c.rs[2]=ratio;}
    s.add_material(c);
    s.add_box(0,D,0,B,0,L,1);
    s.add_fix(-1,1,-1,1,-1e-6,1e-6,true,true,true);
    s.set_gravity(Vec3(0,0,0));
    s.set_damping(20.0);
    s.set_element_type(full);
    s.build();
    std::vector<int> top;
    for(int i=0;i<s.nnode();++i) if(std::fabs(s.ref_pos(i)(2)-L)<1e-9) top.push_back(i);
    const double dt=2e-5, K=2e8*B;
    double ut=0, M=0;
    *Mcr=*My=*uy=*M3=*Mpk=0;
    while(ut<umax)
    {
        ut+=v*dt;
        double u=0; for(int i:top) u+=s.pos(i)(0)-s.ref_pos(i)(0); u/=top.size();
        const double H=K*(ut-u);
        s.clear_loads();
        for(int i:top) s.add_load(i,Vec3(H/top.size(),0,0));
        s.advance(dt);
        M=H/B*L;                                   // base moment per metre width
        if(*Mcr==0 && s.max_damage()>0.05) *Mcr=M;
        if(*My==0 && s.steel_yielded()) {*My=M; *uy=u;}
        if(*M3==0 && u>=0.03*L) *M3=M;
        *Mpk=std::max(*Mpk,M);
    }
    *Mend=M; *rupt=s.bars_ruptured(); *eroded=s.n_eroded();
}

static void test_rc_column()
{
    std::cout<<"RC cantilever 0.5 x 3 m (plane strain), 1 % bars B500B in the outer layer, pushed at the top (10 cm elements, full integration)"<<std::endl;
    // moment-curvature analysis of the section with the same laws (fibre per element column, 10 cm):
    // first yield 571 kNm/m at a curvature of 0.0079 1/m, up to 594 kNm/m beyond; f_t W = 121 kNm/m
    const double My_ref=571e3;
    double Mcr,My,uy,M3,Mpk,Mend; int rupt,er;
    rc_push(0.1,1,0.01,0.1,0.1,&Mcr,&My,&uy,&M3,&Mpk,&Mend,&rupt,&er);
    std::printf("    RC: cracking %.0f kNm/m, first yield %.0f kNm/m at a top drift of %.1f mm, at 3 %% drift %.0f, max %.0f, ruptured %d, eroded %d\n",Mcr/1e3,My/1e3,1e3*uy,M3/1e3,Mpk/1e3,rupt,er);
    check(std::fabs(My/My_ref-1)<0.1,"moment at first yield of the bars within 10 % of the section analysis");
    check(uy>0.017 && uy<0.032,"yield drift: phi_y L^2/3 = 24 mm plus shear (17-32 mm)");
    check(M3>0.95*My && M3<1.25*My,"ductile: at 3 % drift the moment stays between 0.95 and 1.25 M_y");
    double Mcr0,My0,uy0,M30,Mpk0,Mend0; int rupt0,er0;
    rc_push(0.05,0,0.0,0.03,0.05,&Mcr0,&My0,&uy0,&M30,&Mpk0,&Mend0,&rupt0,&er0);
    std::printf("    plain: max %.0f kNm/m (f_t W = 121), end %.0f, eroded %d\n",Mpk0/1e3,Mend0/1e3,er0);
    check(Mpk0>121e3 && Mpk0<0.4*My,"plain concrete: brittle, the peak moment is a fraction of the RC yield moment");
    check(er0>0 && Mend0<0.2*Mpk0,"plain concrete: the section breaks");
}

// ----------------------------------------------------------------------
// threads (OpenMP): the results must not depend on the number of threads
// ----------------------------------------------------------------------

static void threads_run(int nthr,std::vector<double>& state,int* eroded,int* nrigid,double* wd)
{
    // a weak plain concrete wall broken by a rigid block (erosion, fragments, node and
    // ground contact), a reinforced column hit by a second block (bars, unilateral damage)
    std::istringstream in(
        "lattice 0.05 0.05 0.05\n"
        "material 1 concrete 2400 3e9 0.2 2e4 50 2e6 5e3\n"
        "material 2 concrete C30\n"
        "rebar B500B z 1%\n"
        "box 0.5 0.6 0 0.4 0 0.8 1\n"
        "box 1.2 1.4 0 0.2 0 1.0 2\n"
        "fix -1 3 -1 1 -0.001 0.001 xyz\n"
        "ground 0.0\n"
        "contact on\n");
    fem_solid s;
    s.read(in);
    fem_solid::material m=elastic(3,500,1e10,0.3); m.rigid=true; m.kdebris=1e7;
    s.add_material(m);
    s.add_box(0.1,0.3,0.1,0.3,0.3,0.5,3);
    s.add_box(0.85,1.05,0.0,0.2,0.4,0.6,3);
    s.build();
    s.use_threads(nthr);
    for(int k=0;k<s.n_rigid();++k)
    {
        fem_solid::rigid_body& rb=const_cast<fem_solid::rigid_body&>(s.rigid(k));
        rb.V = Vec3(rb.c(0)<0.5 ? 15.0 : 10.0,0,0);
    }
    const int n0=s.n_alive();
    for(int n=0;n<60;++n) s.advance(1e-3);
    state.clear();
    for(int i=0;i<s.nnode();++i) for(int d=0;d<3;++d) {state.push_back(s.pos(i)(d)); state.push_back(s.vel(i)(d));}
    for(int e=0;e<s.nelem();++e) for(int q=0;q<s.ngauss();++q)
    {
        const fem_solid::gpstate& g=s.gp(e,q);
        state.push_back(g.d); state.push_back(g.kt); state.push_back(g.kc);
        for(int k=0;k<3;++k) state.push_back(g.es[k]);
    }
    for(int k=0;k<s.n_rigid();++k) for(int d=0;d<3;++d) {state.push_back(s.rigid(k).c(d)); state.push_back(s.rigid(k).V(d));}
    *eroded=n0-s.n_alive(); *nrigid=s.n_rigid(); *wd=s.dissipated_energy();
    state.push_back(*wd);
}

static void test_threads()
{
    const int T=std::max(2,fem_solid::max_threads());
    std::cout<<"threads: 1 and "<<T<<" threads give bitwise the same results"
             <<(fem_solid::max_threads()>1 ? "" : " (built without OpenMP: one thread)")<<std::endl;
    std::vector<double> a,b;
    int ea,eb,ra,rb; double wa,wb;
    threads_run(1,a,&ea,&ra,&wa);
    threads_run(T,b,&eb,&rb,&wb);
    std::printf("    eroded %d / %d elements, rigid bodies %d / %d, dissipated %.6f / %.6f J, %zu values\n",ea,eb,ra,rb,wa,wb,a.size());
    check(ea>0 && ra>2,"the wall breaks: eroded elements and rigid fragments");
    check(a.size()==b.size() && std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0,"positions, velocities, damage, bar strains, rigid bodies, dissipated energy bitwise equal");
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
    if(w=="all"||w=="snap") test_snap();
    if(w=="all"||w=="patch") test_patch();
    if(w=="all"||w=="presets") test_presets();
    if(w=="all"||w=="settle") test_settle_check();
    if(w=="all"||w=="snapbeam") test_snapped_cantilever();
    if(w=="all"||w=="damping") test_damping();
    if(w=="all"||w=="rigid") test_rigid();
    if(w=="all"||w=="walls") test_walls();
    if(w=="all"||w=="impact") test_impact();
    if(w=="all"||w=="rebar") test_rebar_input();
    if(w=="all"||w=="tie") test_rebar_tension();
    if(w=="all"||w=="rc") test_rc_column();
    if(w=="all"||w=="threads") test_threads();
    std::cout<<(nfail ? "FAILED: " : "all tests passed")<<(nfail? std::to_string(nfail):"")<<std::endl;
    return nfail ? 1 : 0;
}
