// Architect: Hans Bihs
// Standalone verification of the ship module kernels (ship_hull, ship_models): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src ship_test.cpp ../../src/ship_hull.cpp ../../src/ship_models.cpp
//         ../../src/geo_primitive.cpp ../../src/6DOF_actuator_disk.cpp -I../../ThirdParty/eigen-5.0.0
//         -DEIGEN_MPL2_ONLY -o ship_test
// Run:    ./ship_test
#include"ship_hull.h"
#include"ship_models.h"
#include"geo_primitive.h"
#include"6DOF_load.h"
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<cstdio>

static int nfail = 0;
static void check(bool ok, const std::string &what)
{
    std::cout<<(ok ? "  PASS  " : "  FAIL  ")<<what<<std::endl;
    if(!ok) ++nfail;
}

struct mesh
{
    double **x, **y, **z;
    int n=0;
    mesh(int nmax)
    {
        x = new double*[nmax]; y = new double*[nmax]; z = new double*[nmax];
        for(int i=0;i<nmax;++i) {x[i]=new double[3]; y[i]=new double[3]; z[i]=new double[3];}
    }
};

int main()
{
    char s[300];
    
    std::cout<<"triangle clipping"<<std::endl;
    {
        const double a[3]={0,0,0}, b[3]={1,0,0}, c[3]={0,0,1};   // area 0.5 in the x-z plane
        double poly[4][3];
        int np = ship_hull::clip_below(a,b,c,2.0,poly);
        check(np==3 && fabs(ship_hull::polygon_area(poly,np)-0.5)<1e-14, "fully below: whole triangle");
        np = ship_hull::clip_below(a,b,c,-1.0,poly);
        check(ship_hull::polygon_area(poly,np)==0.0, "fully above: nothing");
        np = ship_hull::clip_below(a,b,c,0.5,poly);
        check(np==4 && fabs(ship_hull::polygon_area(poly,np)-0.375)<1e-14, "half height: trapezoid 3/8");
    }
    
    std::cout<<"box hull L=2, B=0.5, H=0.4, CoG at the centre, draft d=0.15"<<std::endl;
    const double L=2.0, B=0.5, H=0.4, d=0.15;
    const double zw = -0.5*H + d;     // still water level in the body frame
    mesh m(100);
    geo_primitive::box(m.x,m.y,m.z,m.n,-0.5*L,0.5*L,-0.5*B,0.5*B,-0.5*H,0.5*H);
    {
        const double S = ship_hull::wetted_surface(m.x,m.y,m.z,m.n,zw);
        const double Sex = L*B + 2.0*(L+B)*d;
        snprintf(s,300,"wetted surface %.12f (exact %.12f)",S,Sex);
        check(fabs(S-Sex)<1e-12, s);
        
        const double S2 = ship_hull::wetted_surface(m.x,m.y,m.z,m.n,zw,true);
        const double S2ex = L*B + 2.0*B*d;
        snprintf(s,300,"2D (no y faces) %.12f (exact %.12f)",S2,S2ex);
        check(fabs(S2-S2ex)<1e-12, s);
        
        double xa,xf;
        ship_hull::waterline_extent(m.x,m.y,m.z,m.n,zw,xa,xf);
        check(fabs(xa+0.5*L)<1e-14 && fabs(xf-0.5*L)<1e-14, "waterline from -L/2 to L/2");
        
        std::vector<double> xs,dx,T;
        ship_hull::draft_strips(m.x,m.y,m.z,m.n,zw,xa,xf,20,xs,dx,T);
        double Tmin=1e9,Tmax=-1e9,sum=0;
        for(size_t i=0;i<T.size();++i){Tmin=std::min(Tmin,T[i]);Tmax=std::max(Tmax,T[i]);sum+=dx[i];}
        snprintf(s,300,"draft strips T = %.12f .. %.12f (exact %.2f), length %.12f",Tmin,Tmax,d,sum);
        check(fabs(Tmin-d)<1e-12 && fabs(Tmax-d)<1e-12 && fabs(sum-L)<1e-12, s);
        
        std::cout<<"cross-flow drag"<<std::endl;
        const double rho=1000.0, Cd=0.8;
        double Y,N;
        ship_models::crossflow(rho,Cd,xs,dx,T,0.3,0.0,Y,N);
        const double Yex = -0.5*rho*Cd*d*L*0.3*0.3;
        snprintf(s,300,"pure sway: Y = %.10f (exact %.10f), N = %.2e",Y,Yex,N);
        check(fabs(Y-Yex)<1e-10 && fabs(N)<1e-10, s);
        
        // pure yaw rate: N = -1/2 rho Cd T r|r| int x^2 |x| dx = -1/2 rho Cd T r|r| (L/2)^4/2
        std::vector<double> xs2,dx2,T2;
        ship_hull::draft_strips(m.x,m.y,m.z,m.n,zw,xa,xf,400,xs2,dx2,T2);
        const double r=0.2;
        ship_models::crossflow(rho,Cd,xs2,dx2,T2,0.0,r,Y,N);
        const double Nex = -0.5*rho*Cd*d*r*fabs(r)*pow(0.5*L,4)/2.0;
        snprintf(s,300,"pure yaw, 400 strips (midpoint rule): N = %.8f (exact %.8f), Y = %.2e",N,Nex,Y);
        check(fabs(N-Nex)<1e-4*fabs(Nex) && fabs(Y)<1e-10, s);
    }
    
    std::cout<<"ITTC-1957 friction"<<std::endl;
    {
        check(fabs(ship_models::cf_ittc57(1.0e7)-0.003)<1e-15, "C_F(1e7) = 0.075/25 = 0.003");
        check(ship_models::cf_ittc57(10.0)==ship_models::cf_ittc57(ship_models::Re_min), "Re below Re_min uses Re_min");
        double Re,CF;
        const double rho=1000.0, nu=1.0e-6, S=1.6, Lw=2.0, k=0.1, u=1.5;
        const double X = ship_models::friction(rho,nu,S,Lw,k,u,Re,CF);
        const double Xex = -0.5*rho*S*1.1*(0.075/pow(log10(3.0e6)-2.0,2))*u*u;
        snprintf(s,300,"X(u=1.5) = %.10f (exact %.10f), Re = %.3e",X,Xex,Re);
        check(fabs(X-Xex)<1e-10 && fabs(Re-3.0e6)<1e-6, s);
        double Xm = ship_models::friction(rho,nu,S,Lw,k,-u,Re,CF);
        check(Xm==-X, "antisymmetric in u");
        check(ship_models::friction(rho,nu,S,Lw,k,0.0,Re,CF)==0.0, "zero at rest");
    }
    
    std::cout<<"propeller open-water model"<<std::endl;
    {
        const double kt[3]={0.5,-0.4,-0.1}, kq[3]={0.06,-0.04,-0.01};
        double J,KT,KQ,Tp,Qp;
        ship_models::propeller(1000.0,10.0,0.2,kt,kq,1.0,J,KT,KQ,Tp,Qp);
        const double Jex=0.5, KTex=0.5-0.2-0.025, KQex=0.06-0.02-0.0025;
        snprintf(s,300,"J = %.6f, KT = %.6f, T = %.6f N (exact %.6f)",J,KT,Tp,1000.0*100.0*pow(0.2,4)*KTex);
        check(fabs(J-Jex)<1e-14 && fabs(KT-KTex)<1e-14 && fabs(KQ-KQex)<1e-14 && fabs(Tp-1000.0*100.0*pow(0.2,4)*KTex)<1e-10
              && fabs(Qp-1000.0*100.0*pow(0.2,5)*KQex)<1e-10, s);
        ship_models::propeller(1000.0,10.0,0.2,kt,kq,0.0,J,KT,KQ,Tp,Qp);
        check(J==0.0 && fabs(Tp-1000.0*100.0*pow(0.2,4)*0.5)<1e-10, "bollard pull J = 0: T = rho n^2 D^4 kt0");
    }
    
    std::cout<<"actuator disk (Hough-Ordway)"<<std::endl;
    {
        sixdof_actuator_disk ad;
        ad.centre = Eigen::Vector3d(1.0,0.5,-0.2);
        ad.axis = Eigen::Vector3d(1.0,0.0,0.0);
        ad.R = 0.1; ad.Rh = 0.02; ad.thickness = 0.04; ad.T = 5.0; ad.Q = 0.3; ad.sense = 1;
        double wa,wt,r; Eigen::Vector3d et;
        check(!ad.weights(Eigen::Vector3d(1.0,0.5,-0.2),wa,wt,et,r), "no weight on the hub");
        check(!ad.weights(Eigen::Vector3d(1.03,0.55,-0.2),wa,wt,et,r), "no weight outside the thickness");
        check(!ad.weights(Eigen::Vector3d(1.0,0.62,-0.2),wa,wt,et,r), "no weight outside the tip");
        // a point on +y (port) of the axis: right-handed rotation about +x moves it to +z
        bool in = ad.weights(Eigen::Vector3d(1.0,0.56,-0.2),wa,wt,et,r);
        check(in && wa>0.0 && wt>0.0 && fabs(r-0.06)<1e-14 && (et-Eigen::Vector3d(0,0,1)).norm()<1e-14, "blade motion: +y -> +z for sense +1");
        // discrete totals on a fine grid, normalised like the coupling: force T along -axis, torque Q
        const int nx=40, ny=120, nz=120;
        const double hx=0.05/nx, hy=0.24/ny, hz=0.24/nz;
        double SA=0, ST=0;
        for(int pass=0; pass<2; ++pass)
        {
            double Fx=0, Mx=0;
            for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) for(int k=0;k<nz;++k)
            {
                const Eigen::Vector3d x(0.975+(i+0.5)*hx, 0.38+(j+0.5)*hy, -0.32+(k+0.5)*hz);
                if(!ad.weights(x,wa,wt,et,r)) continue;
                const double V=hx*hy*hz;
                if(pass==0) {SA+=wa*V; ST+=wt*r*V; continue;}
                const Eigen::Vector3d f = -ad.T*wa/SA*ad.axis + ad.Q*wt/ST*et;
                Fx += f(0)*V;
                Mx += ((x-ad.centre).cross(f))(0)*V;
            }
            if(pass==1)
            {
                snprintf(s,300,"discrete force %.12f (-T = -5), torque about the axis %.12f (Q = 0.3)",Fx,Mx);
                check(fabs(Fx+5.0)<1e-10 && fabs(Mx-0.3)<1e-10, s);
            }
        }
    }
    
    std::cout<<"MMG rudder"<<std::endl;
    {
        ship_models::rudder_param R;
        R.AR=0.01; R.Lambda=1.8; R.xR=-0.5; R.zR=-0.05; R.tR=0.39; R.aH=0.3; R.xH=-0.45; R.eps=1.1; R.kappa=0.5; R.lR=-0.9; R.gammaR=0.4; R.gammaRp=-1.0; R.falpha=0.0;
        double X,Y,N,K,aR,UR,FN;
        ship_models::rudder_mmg(1000.0,R,1.0,0.0,0.0,0.0,0.1,10.0,0.3,0.8,X,Y,N,K,aR,UR,FN);
        check(fabs(Y)<1e-14 && fabs(N)<1e-14 && fabs(X)<1e-14, "straight run, delta = 0: no rudder force");
        const double d=10.0*3.14159265358979/180.0;
        ship_models::rudder_mmg(1000.0,R,1.0,0.0,0.0,d,0.1,10.0,0.3,0.8,X,Y,N,K,aR,UR,FN);
        snprintf(s,300,"delta = +10 deg: Y = %.5f N > 0, N = %.5f Nm < 0 (turn to starboard), X = %.5f N < 0",Y,N,X);
        check(Y>0.0 && N<0.0 && X<0.0 && fabs(aR-d)<1e-14, s);
        double X2,Y2,N2,K2;
        ship_models::rudder_mmg(1000.0,R,1.0,0.0,0.0,-d,0.1,10.0,0.3,0.8,X2,Y2,N2,K2,aR,UR,FN);
        check(fabs(Y2+Y)<1e-12 && fabs(N2+N)<1e-12 && fabs(X2-X)<1e-12, "antisymmetric in delta");
        // propeller slipstream accelerates the rudder inflow
        double UR0;
        ship_models::rudder_mmg(1000.0,R,1.0,0.0,0.0,d,0.0,0.0,0.0,0.8,X,Y,N,K,aR,UR0,FN);
        check(UR>UR0 && fabs(UR0-1.1*0.8)<1e-14, "slipstream: U_R with propeller > without (= eps uP)");
        // ship drifting to port (v > 0) at delta = 0: flow on the rudder from port, force to starboard
        ship_models::rudder_mmg(1000.0,R,1.0,0.1,0.0,0.0,0.1,10.0,0.3,0.8,X,Y,N,K,aR,UR,FN);
        check(Y<0.0, "drift to port at delta = 0: rudder side force to starboard (course stability)");
    }
    
    std::cout<<"MMG rudder options"<<std::endl;
    {
        ship_models::rudder_param R;
        R.AR=0.0539; R.Lambda=0.345*0.345/0.0539; R.xR=-3.75; R.zR=0.0; R.tR=0.387; R.aH=0.312; R.xH=-3.248; R.eps=1.09; R.kappa=0.5;
        R.lR=-0.71*7.0; R.gammaR=0.395; R.gammaRp=0.640; R.falpha=2.747;
        double X,Y,N,K,aR,UR,FN,UR2,FN2;
        // falpha given: F_N = 1/2 rho A_R U_R^2 falpha sin(alpha_R)
        ship_models::rudder_mmg(1000.0,R,1.0,0.0,0.0,0.2,0.0,0.0,0.0,0.6,X,Y,N,K,aR,UR,FN);
        check(fabs(FN - 0.5*1000.0*0.0539*UR*UR*2.747*sin(0.2))<1e-10, "rudder_falpha overrides 6.13 Lambda/(Lambda+2.25)");
        // asymmetric flow straightening: beta_R > 0 (drift to port in the ship frame) uses gammaRp
        ship_models::rudder_mmg(1000.0,R,1.0, 0.05,0.0,0.0,0.0,0.0,0.0,0.6,X,Y,N,K,aR,UR,FN);
        ship_models::rudder_mmg(1000.0,R,1.0,-0.05,0.0,0.0,0.0,0.0,0.0,0.6,X,Y,N,K,aR,UR2,FN2);
        const double U = sqrt(1.0+0.0025), b = atan2(0.05,1.0);
        check(fabs(UR - sqrt(1.09*0.6*1.09*0.6 + pow(U*0.640*b,2)))<1e-12 && fabs(UR2 - sqrt(1.09*0.6*1.09*0.6 + pow(U*0.395*b,2)))<1e-12,
              "rudder_gamma: gamma_R 0.640 for beta_R > 0, 0.395 for beta_R < 0 (KVLCC2)");
    }
    
    std::cout<<"MMG hull (Yasukawa & Yoshimura 2015, KVLCC2)"<<std::endl;
    {
        const double c[17] = {0.022,-0.040,0.002,0.011,0.771,-0.315,0.083,-1.607,0.379,-0.391,0.008,-0.137,-0.049,-0.030,-0.294,0.055,-0.013};
        const double rho=1000.0, L=7.0, d=0.46;
        double X,Y,N;
        // straight run: only the resistance -R0
        ship_models::mmg_hull(rho,L,d,c,1.2,0.0,0.0,0.01,X,Y,N);
        check(fabs(X + 0.5*rho*L*d*1.44*0.022)<1e-10 && Y==0.0 && N==0.0, "straight run: X = -1/2 rho L d U^2 R0', Y = N = 0");
        // general state against the nondimensional form
        const double u=1.1, vm=-0.12, r=0.03, U=sqrt(u*u+vm*vm), vp=vm/U, rp=r*L/U;
        ship_models::mmg_hull(rho,L,d,c,u,vm,r,0.01,X,Y,N);
        const double q=0.5*rho*L*d*U*U;
        const double Xe=q*(-c[0]+c[1]*vp*vp+c[2]*vp*rp+c[3]*rp*rp+c[4]*pow(vp,4));
        const double Ye=q*(c[5]*vp+c[6]*rp+c[7]*pow(vp,3)+c[8]*vp*vp*rp+c[9]*vp*rp*rp+c[10]*pow(rp,3));
        const double Ne=q*L*(c[11]*vp+c[12]*rp+c[13]*pow(vp,3)+c[14]*vp*vp*rp+c[15]*vp*rp*rp+c[16]*pow(rp,3));
        snprintf(s,300,"polynomial form = nondimensional form: X %.6f/%.6f, Y %.6f/%.6f, N %.6f/%.6f",X,Xe,Y,Ye,N,Ne);
        check(fabs(X-Xe)<1e-10 && fabs(Y-Ye)<1e-10 && fabs(N-Ne)<1e-10, s);
        // frame independence: (v, r) -> (-v, -r) flips Y and N, keeps X
        double X2,Y2,N2;
        ship_models::mmg_hull(rho,L,d,c,u,-vm,-r,0.01,X2,Y2,N2);
        check(fabs(X2-X)<1e-12 && fabs(Y2+Y)<1e-12 && fabs(N2+N)<1e-12, "ship frame = MMG frame: Y, N odd and X even in (v, r)");
        // course stability terms: sway to port (vm > 0) gives a force to starboard, yaw damping
        ship_models::mmg_hull(rho,L,d,c,1.0,0.05,0.0,0.01,X,Y,N);
        check(Y<0.0, "Y'v < 0: lateral damping");
        ship_models::mmg_hull(rho,L,d,c,1.0,0.0,0.02,0.01,X,Y,N);
        check(N<0.0, "N'r < 0: yaw damping");
        // finite at rest with a yaw rate (U limited to Umin in the denominators)
        ship_models::mmg_hull(rho,L,d,c,0.0,0.0,0.01,0.05,X,Y,N);
        check(std::isfinite(X) && std::isfinite(Y) && std::isfinite(N), "finite at U = 0");
    }
    
    std::cout<<"MMG wake"<<std::endl;
    {
        check(fabs(ship_models::mmg_wake(0.40,2.0,1.6,1.1,0.0)-0.40)<1e-15, "beta_P = 0: wP = wP0");
        const double w1 = ship_models::mmg_wake(0.40,2.0,1.6,1.1,0.3), w2 = ship_models::mmg_wake(0.40,2.0,1.6,1.1,-0.3);
        const double e1 = 1.0-0.6*(1.0+(1.0-exp(-0.6))*0.6), e2 = 1.0-0.6*(1.0+(1.0-exp(-0.6))*0.1);
        check(fabs(w1-e1)<1e-14 && fabs(w2-e2)<1e-14, "C2 = 1.6 for beta_P > 0, 1.1 for beta_P < 0");
    }
    
    std::cout<<"roll damping"<<std::endl;
    check(ship_models::roll_damping(10.0,5.0,-0.2)==10.0*0.2+5.0*0.04, "K = -B44 p - B44q |p| p");
    
    std::cout<<(nfail ? "FAILED: " : "all tests passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
