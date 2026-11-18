#include <algorithm>
// Standalone check of the rheology core (src/rheology_core.{h,cpp}): no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src rheology_core_test.cpp ../../src/rheology_core.cpp -o rheology_core_test
// Run:    ./rheology_core_test
// Checks: Herschel-Bulkley plug flow against the Bingham closed form and its inverse, the Newtonian and
// power-law limits, the backward-Euler friction step (Newtonian exact, Coulomb stop and exact speed, pore
// pressure, residuals of Voellmy, Pouliquen-Forterre, mu(I), HB, HBC), the partially regularised mu(I)
// (neutral points, C <= 0 on the regularised branch and on [I_N, I_2N], continuity), the Papanastasiou
// viscosity, the HBC yield stress without tension and the pore pressure relaxation, mu(J) / mu(K) (Boyer et al. 2011,
// Trulsson et al. 2012) and the dilatancy friction tan(psi).
// Architect: Hans Bihs
#include"rheology_core.h"
#include <cstdio>
#include <cmath>
using namespace rheo;
static int fails=0;
#define CHECK(c,msg) do{ if(!(c)){ printf("FAIL: %s\n",msg); ++fails;} else printf("ok:   %s\n",msg);}while(0)
int main(){
  // 1. Bingham plug flow vs closed form U = tau h/(6 mu)(2 - 3r + r^3)
  { double ty=50, mu=10, h=0.05, tb=120; double r=ty/tb;
    double Ua=tb*h/(6*mu)*(2-3*r+r*r*r), Un=hb_mean_speed(tb,ty,mu,1.0,h);
    CHECK(fabs(Un-Ua)<1e-12*Ua,"HB n=1 equals the Bingham closed form");
    double t2=hb_basal_stress(Un,ty,mu,1.0,h); CHECK(fabs(t2-tb)<1e-9*tb,"HB stress inverts the mean speed");}
  // 2. Newtonian: tau = 3 mu U/h
  { double mu=2, h=0.1, U=0.3; CHECK(fabs(hb_basal_stress(U,0,mu,1.0,h)-3*mu*U/h)<1e-9,"tau_y=0, n=1: 3 mu U/h");}
  // 3. power law n=0.4: U = a n h^(1+1/n)/(2n+1), a=(tau/(K h))^(1/n)
  { double K=20,n=0.4,h=0.05,tb=30; double a=pow(tb/(K*h),1/n); double Ua=a*n*pow(h,1+1/n)/(2*n+1);
    CHECK(fabs(hb_mean_speed(tb,0,K,n,h)-Ua)<1e-12*Ua,"power-law mean speed");}
  // 4. implicit Newtonian = analytic backward Euler
  { basal_params q; q.law=BASAL_HB; q.tauy=0; q.K=2; q.n=1; basal_state s; s.h=0.1; s.rho=1000; s.gn=9.81; s.sigma_n=s.rho*s.gn*s.h;
    double dt=0.01, Us=1.0; double c=dt/(s.rho*s.h); double ex=Us/(1+c*3*q.K/s.h);
    CHECK(fabs(implicit_speed(q,s,Us,dt)-ex)<1e-10,"implicit Newtonian = backward Euler");}
  // 5. Coulomb: stops below the threshold, exact otherwise
  { basal_params q; q.law=BASAL_COULOMB; q.tandelta=tan(deg2rad(30)); basal_state s; s.h=0.1; s.rho=2000; s.gn=9.81; s.sigma_n=s.rho*s.gn*s.h;
    double dt=0.01, c=dt/(s.rho*s.h), ts=q.tandelta*s.sigma_n;
    CHECK(implicit_speed(q,s,0.9*c*ts,dt)==0.0,"Coulomb: stop below threshold");
    CHECK(fabs(implicit_speed(q,s,2*c*ts,dt)-c*ts)<1e-12,"Coulomb: s = U* - c tau");
    s.lambda=0.5; CHECK(fabs(basal_stress(q,s,1.0)-0.5*ts)<1e-9,"pore pressure halves the friction");}
  // 6. Voellmy closed form residual
  { basal_params q; q.law=BASAL_VOELLMY; q.tandelta=0.2; q.xi=500; basal_state s; s.h=1; s.rho=1500; s.gn=9.81; s.sigma_n=s.rho*s.gn*s.h;
    double dt=0.05, Us=20, sp=implicit_speed(q,s,Us,dt); double res=sp+dt/(s.rho*s.h)*basal_stress(q,s,sp)-Us;
    CHECK(fabs(res)<1e-9,"Voellmy implicit residual");}
  // 7. Pouliquen-Forterre and mu(I) implicit residuals, HB residual
  { basal_params q; q.law=BASAL_PF; q.tz1=tan(deg2rad(21)); q.tz2=tan(deg2rad(30.7)); q.tz3=tan(deg2rad(22.2)); q.beta=0.136; q.L=6.5e-4;
    basal_state s; s.h=0.005; s.rho=1500; s.gn=9.81*cos(deg2rad(28)); s.sigma_n=s.rho*s.gn*s.h;
    double dt=1e-3, Us=0.5, sp=implicit_speed(q,s,Us,dt); double res=sp+dt/(s.rho*s.h)*basal_stress(q,s,sp)-Us;
    CHECK(fabs(res)<1e-9,"Pouliquen-Forterre implicit residual");
    q.law=BASAL_MUI; sp=implicit_speed(q,s,Us,dt); res=sp+dt/(s.rho*s.h)*basal_stress(q,s,sp)-Us; CHECK(fabs(res)<1e-9,"mu(I) implicit residual");
    q.law=BASAL_HB; q.tauy=20; q.K=5; q.n=0.4; sp=implicit_speed(q,s,Us,dt); res=sp+dt/(s.rho*s.h)*basal_stress(q,s,sp)-Us; CHECK(fabs(res)<1e-7,"HB implicit residual");
    q.law=BASAL_HBC; q.tandelta=0.3; sp=implicit_speed(q,s,Us,dt); res=sp+dt/(s.rho*s.h)*basal_stress(q,s,sp)-Us; CHECK(fabs(res)<1e-7,"HBC implicit residual");}
  // 7b. sub-stepped friction of a long step (semi-implicit SFLOW): du/dt = a - tau(u)/(rho h) against a fine RK4 solution
  { basal_params q; q.law=BASAL_MUI; basal_state s; s.h=0.005; s.rho=1500; s.gn=9.81*cos(deg2rad(25)); s.sigma_n=s.rho*s.gn*s.h;
    const double a=9.81*sin(deg2rad(25)), dt=0.04;
    auto f=[&](double u){ return a - basal_stress(q,s,std::max(u,0.0))/(s.rho*s.h); };
    double u=0, ue=0, ub=0, emax=0, ebmax=0;
    for(int k=0;k<25;++k){ u=implicit_speed_sub(q,s,u,u+dt*a,dt,8); ub=implicit_speed(q,s,ub+dt*a,dt);
      for(int m=0;m<4000;++m){ double h=dt/4000, k1=f(ue), k2=f(ue+0.5*h*k1), k3=f(ue+0.5*h*k2), k4=f(ue+h*k3); ue+=h*(k1+2*k2+2*k3+k4)/6; }
      emax=fmax(emax,fabs(u-ue)/0.0869); ebmax=fmax(ebmax,fabs(ub-ue)/0.0869); }
    printf("      mu(I) layer, dt 0.04 s: max error sub-stepped %.3f, single step %.3f (of the steady speed)\n",emax,ebmax);
    CHECK(emax<0.02 && emax<0.25*ebmax,"sub-stepped friction follows the acceleration");
    q.law=BASAL_COULOMB; q.tandelta=tan(deg2rad(30)); s.sigma_n=s.rho*9.81*s.h; s.gn=9.81;
    CHECK(implicit_speed_sub(q,s,1.0,0.5,dt,8)==implicit_speed(q,s,0.5,dt),"Coulomb: sub-stepping = the exact single step");}
  // 8. steady uniform PF flow: h_stop law consistency, mu_start(h) > mu_dyn(h, beta)
  // 9. mu(I) partial regularisation (glass beads, Jop et al.: mus=tan20.9, mu2=tan32.76, I0=0.279)
  { law3D_params q; q.law=LAW_MUI; q.mus=tan(deg2rad(20.9)); q.mu2=tan(deg2rad(32.76)); q.I0=0.279; q.alpha=1.9;
    double I2=mu_I_setup(q); printf("      I_N = %.4e, I_2N = %.4e, A = %.4e\n",q.I_N,I2,q.A);
    CHECK(q.I_N>0 && I2>q.I_N,"neutral points found");
    double mx=-1e9; for(int k=0;k<4000;++k){ double I=pow(10.0,-12+k*(log10(q.I_N)+12)/4000.0)*0.999999; double mu=mu_I(q,I);
       double h=I*1e-6; double nu=I*(mu_I(q,I+h)-mu_I(q,I-h))/(2*h)/mu; mx=fmax(mx,well_posedness(mu,nu));}
    printf("      max C on the regularised branch = %.3e\n",mx);
    CHECK(mx<=1e-9,"regularised branch well posed (C <= 0)");
    CHECK(fabs(mu_I(q,q.I_N*(1-1e-12))-mu_I_plain(q,q.I_N))<1e-9,"continuous at I_N");
    double mxp=-1e9; for(int k=0;k<1000;++k){ double I=q.I_N*pow(I2/q.I_N,(k+0.5)/1000.0); mxp=fmax(mxp,well_posedness(mu_I_plain(q,I),I*(q.mu2-q.mus)*q.I0/pow(q.I0+I,2)/mu_I_plain(q,I)));}
    CHECK(mxp<=1e-12,"plain law well posed on [I_N, I_2N]");
    // viscosity: finite, capped, consistent with tau = mu(I) p
    double p=1000, g=10; double mu=viscosity(q,g,p); double I=inertial_number(g,p,q.d,q.rho_s);
    CHECK(fabs(mu*g - mu_I(q,I)*p)<1e-9*p,"mu_eff gamma = mu(I) p'");}
  // 10. HB Papanastasiou: tau -> tau_y as gamma -> 0
  { law3D_params q; q.law=LAW_HB; q.tauy=40; q.K=10; q.n=1; q.m=1e4; q.mu_max=1e9;
    double g=1e-2; CHECK(fabs(viscosity(q,g,0)*g-(40*(1-exp(-1e4*g))+10*g))<1e-9,"HB Papanastasiou stress");
    q.law=LAW_HBC; q.tanphi=0.5; CHECK(fabs(yield_stress(q,-50)-40)<1e-12 && fabs(yield_stress(q,100)-90)<1e-12,"HBC yield stress, no tension");}
  // 11. pore pressure relaxation
  { double l=pore_relax(0.8,0.4,1e-3,0.1,1.0); double ex=0.4+0.4*exp(-rheo::PI_*rheo::PI_*1e-3/(4*0.01)); CHECK(fabs(l-ex)<1e-14,"pore pressure relaxation");}
  // 12. mu(J) (Boyer, Guazzelli & Pouliquen 2011: mu1 0.32, mu2 0.7, J0 0.005, phi_c 0.585)
  { law3D_params q; q.law=LAW_MUJ; q.mus=0.32; q.mu2=0.7; q.J0=0.005; q.phic=0.585; q.eta_f=2.15; q.d=1.1e-3; q.rho_s=1050; q.mu_max=1e12;
    double J=0.01; double ex=0.32+(0.7-0.32)/(1+0.005/J)+J+2.5*0.585*sqrt(J);
    CHECK(fabs(mu_J(q,J,0.0)-ex)<1e-14,"mu(J) = Boyer et al. (2011)");
    CHECK(fabs(mu_J(q,1e-14,0.0)-0.32)<1e-6,"mu(J -> 0) = mu_s");
    // viscosity: tau = mu(J) p' at the J of (gamma, p')
    double p=200, g=3; double mu=viscosity(q,g,p); double Jg=viscous_number(g,p,q.eta_f);
    CHECK(fabs(mu*g - mu_J(q,Jg,inertial_number(g,p,q.d,q.rho_s))*p)<1e-9*p,"mu_eff gamma = mu(J) p'");
    // p' = 0: the fluid viscosity
    CHECK(fabs(viscosity(q,g,0.0)-q.eta_f)<1e-12,"mu(J) at p' = 0: eta_f");
    // large J: stress ratio -> J, mu_eff -> eta_f (the suspension limit of the fluid term)
    CHECK(fabs(viscosity(q,g,1e-9)-q.eta_f)/q.eta_f<1e-3,"mu(J) at J >> 1: mu_eff -> eta_f");
    // mixed number: aK = 0 is J, aK > 0 raises the contact friction in the inertial regime
    double I=inertial_number(g,p,q.d,q.rho_s); q.aK=0.635;
    double K=Jg+0.635*I*I; CHECK(mu_J(q,Jg,I)>mu_J(q,Jg,0.0) && fabs(mu_J(q,Jg,I)-(0.32+0.38/(1+0.005/K)+Jg+2.5*0.585*sqrt(Jg)))<1e-14,"mu(K): K = J + aK I^2 in the contact part");
    q.aK=0.0;
    CHECK(fabs(yield_stress(q,100.0)-32.0)<1e-12,"mu(J) yield stress mu_s p'");
    // dilatancy: tan(psi) = K3 (phi - phi_c/(1 + J^1/2)), added to the friction, never below zero friction
    CHECK(fabs(phi_eq(0.585,0.0)-0.585)<1e-15 && fabs(tan_psi(4.09,0.6,0.585,0.0)-4.09*0.015)<1e-14,"phi_eq(0) = phi_c, tan(psi)");
    double tp=0.1; double mud=viscosity(q,g,p,tp);
    CHECK(fabs(mud*g-(mu_J(q,Jg,I)+tp)*p)<1e-9*p,"dilatancy friction mu(J) + tan(psi)");
    CHECK(viscosity(q,g,p,-10.0)==q.mu_min,"contracting: friction floored at 0");
    CHECK(viscosity(q,g,p,0.0)==viscosity(q,g,p),"tan(psi) = 0: the plain law (bitwise)");
    law3D_params qi; qi.law=LAW_MUI; qi.mus=0.38; qi.mu2=0.64; qi.alpha=0.0; double mi=viscosity(qi,g,p,tp);
    CHECK(fabs(mi*g-(mu_I(qi,inertial_number(g,p,qi.d,qi.rho_s))+tp)*p)<1e-9*p,"dilatancy friction mu(I) + tan(psi)");}
  printf(fails?"\n%d FAILED\n":"\nall passed\n",fails); return fails;
}
