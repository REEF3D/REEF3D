// Standalone verification of the ship module kernels (ship_hull, ship_models): no MPI, no REEF3D.
// Build:  g++ -O2 -std=c++20 -I../../src ship_test.cpp ../../src/ship_hull.cpp ../../src/ship_models.cpp
//         ../../src/geo_primitive.cpp -o ship_test
// Run:    ./ship_test
#include"ship_hull.h"
#include"ship_models.h"
#include"geo_primitive.h"
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
    
    std::cout<<"roll damping"<<std::endl;
    check(ship_models::roll_damping(10.0,5.0,-0.2)==10.0*0.2+5.0*0.04, "K = -B44 p - B44q |p| p");
    
    std::cout<<(nfail ? "FAILED: " : "all tests passed")<<(nfail ? std::to_string(nfail) : "")<<std::endl;
    return nfail ? 1 : 0;
}
