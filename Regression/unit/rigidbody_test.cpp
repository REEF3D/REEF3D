// Standalone verification of the 6DOF rigid-body core (sixdof_rigidbody): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src
//         rigidbody_test.cpp ../../src/6DOF_rigidbody.cpp -o rigidbody_test
// Run:    ./rigidbody_test
//
// Tests: quaternion/Euler round trip, constant force (exact for all schemes), spring-mass
// oscillator (period and convergence order), torque-free asymmetric top (conservation of
// angular momentum and energy, convergence order), steady spin, DOF modes and damping.
#include"6DOF_rigidbody.h"
#include<iostream>
#include<cmath>
#include<functional>
#include<string>
#include<cstdio>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

enum scheme {RK2, RK3, RKLS3, RK4};
static const char *name(scheme s) {return s==RK2?"RK2":s==RK3?"RK3":s==RKLS3?"RKLS3":"RK4";}
static int stages(scheme s) {return s==RK2 ? 2 : s==RK4 ? 4 : 3;}

// external load as a function of the body state: X, Y, Z, K, M, N
typedef std::function<void(const sixdof_rigidbody&, double*)> loadfunc;

static void step(sixdof_rigidbody &b, scheme s, double dt, const loadfunc &load)
{
    static const double gam[3] = {8.0/15.0, 5.0/12.0, 3.0/4.0};
    static const double zet[3] = {0.0, -17.0/60.0, -5.0/12.0};

    for(int iter=0; iter<stages(s); ++iter)
    {
        double Fext[6];
        load(b,Fext);
        b.assemble_loads(Fext);
        b.derivatives_trans();
        b.derivatives_rot();

        if(s==RK2)   b.stage_rk2(iter,dt);
        if(s==RK3)   b.stage_rk3(iter,dt);
        if(s==RKLS3) b.stage_rkls3(iter,gam[iter],zet[iter],dt);
        if(s==RK4)   b.stage_rk4(iter,dt);

        b.quat_matrices();
        b.update_omega();
    }
    b.save_history(dt);
}

static void setup(sixdof_rigidbody &b, double m, const Eigen::Matrix3d &I)
{
    b.reset();
    b.mass = m;
    b.I = I;
    b.quaternion_from_euler();
    b.quat_matrices();
    b.init_history();
    b.update_omega();
}

static void noload(const sixdof_rigidbody&, double *F) {for(int n=0;n<6;++n) F[n]=0.0;}

// ---------------------------------------------------------------------------------------------

static void test_euler()
{
    std::cout<<"quaternion / Euler angles"<<std::endl;
    sixdof_rigidbody b;
    double maxerr=0.0, orth=0.0;
    const double ang[4][3] = {{0.3,-0.2,1.1},{-1.0,0.7,-2.5},{0.0,0.0,0.0},{2.0,-1.2,3.0}};
    for(int k=0; k<4; ++k)
    {
        b.phi=ang[k][0]; b.theta=ang[k][1]; b.psi=ang[k][2];
        b.quaternion_from_euler();
        b.quat_matrices();
        orth = std::max(orth, (b.R*b.R.transpose() - Eigen::Matrix3d::Identity()).norm());
        orth = std::max(orth, std::fabs(b.R.determinant()-1.0));
        b.euler_angles();
        maxerr = std::max(maxerr, std::fabs(b.phi-ang[k][0]) + std::fabs(b.theta-ang[k][1]) + std::fabs(b.psi-ang[k][2]));
    }
    char s[200];
    snprintf(s,200,"round trip error %.2e, |R R^T - 1| + |det R - 1| = %.2e",maxerr,orth);
    check(maxerr<1e-12 && orth<1e-12, s);

    // R maps the body x-axis for a pure yaw psi to (cos psi, sin psi, 0)
    b.phi=0.0; b.theta=0.0; b.psi=0.5;
    b.quaternion_from_euler(); b.quat_matrices();
    Eigen::Vector3d x = b.R*Eigen::Vector3d(1,0,0);
    check(std::fabs(x(0)-cos(0.5))<1e-12 && std::fabs(x(1)-sin(0.5))<1e-12, "yaw rotates the body x-axis counter-clockwise");
}

static void test_constant_force()
{
    std::cout<<"constant force (exact for all schemes)"<<std::endl;
    const double m=2.0, dt=0.01, T=1.0;
    const Eigen::Vector3d F(1.0,-2.0,-19.62), v0(0.5,0.0,1.0), c0(1.0,2.0,3.0);
    for(scheme s : {RK2,RK3,RKLS3,RK4})
    {
        sixdof_rigidbody b;
        setup(b,m,Eigen::Matrix3d::Identity());
        b.c = c0; b.p = m*v0; b.init_history();
        loadfunc load = [&](const sixdof_rigidbody&, double *Fe){Fe[0]=F(0);Fe[1]=F(1);Fe[2]=F(2);Fe[3]=Fe[4]=Fe[5]=0.0;};
        const int n = int(T/dt+0.5);
        for(int k=0;k<n;++k) step(b,s,dt,load);
        Eigen::Vector3d cex = c0 + v0*T + 0.5*F/m*T*T;
        double err = (b.c-cex).norm() + (b.p/m - (v0+F/m*T)).norm();
        char str[200]; snprintf(str,200,"%-5s error %.2e",name(s),err);
        check(err<1e-10, str);
    }
}

// spring-mass oscillator in heave: error at t=T for dt and dt/2
static double oscillator(scheme s, double dt, double *period=nullptr, double T=3.0)
{
    const double m=10.0, k=40.0, z0=0.1;
    sixdof_rigidbody b;
    setup(b,m,Eigen::Matrix3d::Identity());
    b.c(2) = z0; b.init_history();
    loadfunc load = [&](const sixdof_rigidbody &r, double *Fe){Fe[0]=Fe[1]=Fe[3]=Fe[4]=Fe[5]=0.0; Fe[2] = -k*r.c(2);};
    const int n = int(T/dt+0.5);
    double zold=b.c(2), t=0.0, t1=-1.0, t2=-1.0;
    for(int i=0;i<n;++i)
    {
        step(b,s,dt,load); t+=dt;
        if(period && zold<0.0 && b.c(2)>=0.0)
        {
            double tc = t - dt*b.c(2)/(b.c(2)-zold);
            if(t1<0.0) t1=tc; else if(t2<0.0) t2=tc;
        }
        zold=b.c(2);
    }
    if(period) *period = t2-t1;
    const double w = sqrt(k/m);
    return std::fabs(b.c(2) - z0*cos(w*T));
}

static void test_oscillator()
{
    std::cout<<"spring-mass oscillator"<<std::endl;
    double Tp;
    oscillator(RK3,1.0e-3,&Tp,8.0);
    char str[200];
    snprintf(str,200,"RK3 period %.6f s (exact %.6f)",Tp,2.0*M_PI/2.0);
    check(std::fabs(Tp-M_PI)<1e-4, str);

    const double order_min[4] = {1.9,2.9,2.9,3.9};
    for(scheme s : {RK2,RK3,RKLS3,RK4})
    {
        double e1 = oscillator(s,0.02), e2 = oscillator(s,0.01);
        double ord = log(e1/e2)/log(2.0);
        snprintf(str,200,"%-5s order %.2f (errors %.2e, %.2e)",name(s),ord,e1,e2);
        check(ord>order_min[s], str);
    }
}

// torque-free asymmetric top (intermediate axis): |h| and the kinetic energy are invariants
static void top(scheme s, double dt, double T, Eigen::Vector3d &hI, double &dh, double &dE)
{
    Eigen::Matrix3d I = Eigen::Vector3d(1.0,2.0,3.0).asDiagonal();
    sixdof_rigidbody b;
    setup(b,1.0,I);
    b.h = I*Eigen::Vector3d(0.01,2.0,0.01);
    b.init_history(); b.update_omega();
    const double h0 = b.h.norm(), E0 = 0.5*b.h.dot(I.inverse()*b.h);
    const int n = int(T/dt+0.5);
    dh=0.0; dE=0.0;
    for(int k=0;k<n;++k)
    {
        step(b,s,dt,noload);
        dh = std::max(dh, std::fabs(b.h.norm()-h0)/h0);
        dE = std::max(dE, std::fabs(0.5*b.h.dot(I.inverse()*b.h)-E0)/E0);
    }
    hI = b.R*b.h;   // angular momentum in the inertial frame: constant
}

static void test_top()
{
    std::cout<<"torque-free asymmetric top"<<std::endl;
    Eigen::Vector3d h1, h2;
    double dh, dE;
    char str[200];
    // RK3 and RKLS3 combine the unnormalised quaternion stage values and normalise only for the
    // stage kinematics; normalising the combined stage values would reduce them to second order.
    const double order_min[4] = {1.9,2.9,2.9,3.9};
    for(scheme s : {RK2,RK3,RKLS3,RK4})
    {
        top(s,4.0e-3,4.0,h1,dh,dE);
        top(s,2.0e-3,4.0,h2,dh,dE);
        Eigen::Vector3d hI0 = Eigen::Vector3d(1.0,2.0,3.0).asDiagonal()*Eigen::Vector3d(0.01,2.0,0.01);
        double e1=(h1-hI0).norm(), e2=(h2-hI0).norm();
        double ord = log(e1/e2)/log(2.0);
        snprintf(str,200,"%-5s inertial h drift %.2e -> %.2e (order %.2f), max rel. |h| error %.1e, energy %.1e",name(s),e1,e2,ord,dh,dE);
        check(ord>order_min[s] && dh<1e-4 && dE<1e-4, str);
    }
}

static void test_spin()
{
    std::cout<<"steady spin about the z-axis (RK4)"<<std::endl;
    Eigen::Matrix3d I = Eigen::Vector3d(1.0,2.0,3.0).asDiagonal();
    sixdof_rigidbody b;
    setup(b,1.0,I);
    const double w=0.7, dt=0.01, T=2.0;
    b.h = I*Eigen::Vector3d(0,0,w); b.init_history(); b.update_omega();
    for(int k=0;k<int(T/dt+0.5);++k) step(b,RK4,dt,noload);
    b.euler_angles();
    char str[200]; snprintf(str,200,"psi %.10f (exact %.10f), omega_I z %.10f",b.psi,w*T,b.omega_I(2));
    check(std::fabs(b.psi-w*T)<1e-10 && std::fabs(b.omega_I(2)-w)<1e-12, str);
}

static void test_dofs()
{
    std::cout<<"DOF modes and damping"<<std::endl;
    sixdof_rigidbody b;
    setup(b,2.0,Eigen::Matrix3d::Identity());
    b.dof[0]=1; b.dof[1]=0; b.dof[2]=2; b.dof[3]=0; b.dof[4]=1; b.dof[5]=0;
    b.Cdamp_t[0]=3.0; b.Cdamp_r[1]=5.0;
    b.p << 4.0, 1.0, 1.0;
    b.omega_I << 0.0, 0.2, 0.0;
    const double Fe[6] = {10.0, 10.0, 10.0, 1.0, 1.0, 1.0};
    b.assemble_loads(Fe);
    check(std::fabs(b.F(0)-(10.0-3.0*4.0/2.0))<1e-14 && b.F(1)==0.0 && b.F(2)==0.0, "force: free with damping, fixed and prescribed DOFs get zero");
    check(b.M(0)==0.0 && std::fabs(b.M(1)-(1.0-5.0*0.2))<1e-14 && b.M(2)==0.0, "moment: free with damping, fixed DOFs get zero");
    b.dc << 0.0, 0.0, 0.3;
    Eigen::Matrix<double,6,1> u;
    b.velocity(u);
    check(u(0)==2.0 && u(1)==0.0 && u(2)==0.3, "velocity: momentum for free, motionext for prescribed DOFs");
    b.twoD=true;
    check(b.fixed(1) && b.fixed(3) && b.fixed(5) && !b.fixed(0) && b.fixed(2) && !b.fixed(4), "fixed(): 2D and DOF modes");
}

int main()
{
    test_euler();
    test_constant_force();
    test_oscillator();
    test_top();
    test_spin();
    test_dofs();
    std::cout<<(nfail ? "FAILED: " : "all tests passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
