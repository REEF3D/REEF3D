// Architect: Hans Bihs
// Reference for 05_uniform_layer_laws: speed of an infinite uniform layer on an incline,
// du/dt = g_x - tau_b(u)/(rho h), with the basal closure of the rheology core, integrated with
// backward-Euler steps of 1e-6 s (rheo::implicit_speed). Prints t u every 0.05 s.
// Build: g++ -O2 -std=c++17 -I../../../src uniform_ode.cpp ../../../src/rheology_core.cpp -o uniform_ode
// Run:   ./uniform_ode law theta h rho tend [law parameters as in ctrl.txt: see rheo_eval.py]
#include"rheology_core.h"
#include<cstdio>
#include<cstdlib>
#include<cmath>
using namespace rheo;
int main(int argc, char **argv)
{
    basal_params q; q.law = atoi(argv[1]);
    const double th = deg2rad(atof(argv[2])), h = atof(argv[3]), rho = atof(argv[4]), tend = atof(argv[5]);
    double *a = new double[16]; for(int k=0;k<16;++k) a[k] = (argc>6+k) ? atof(argv[6+k]) : 0.0;
    if(q.law==BASAL_VOELLMY) { q.tandelta = tan(deg2rad(a[0])); q.xi = a[1]; }
    if(q.law==BASAL_PF)  { q.tz1=tan(deg2rad(a[0])); q.tz2=tan(deg2rad(a[1])); q.tz3=tan(deg2rad(a[2])); q.beta=a[3]; q.L=a[4]; }
    if(q.law==BASAL_MUI) { q.mus=tan(deg2rad(a[0])); q.mu2=tan(deg2rad(a[1])); q.I0=a[2]; q.d=a[3]; q.phis=a[4]; }
    if(q.law==BASAL_HB || q.law==BASAL_HBC) { q.tauy=a[0]; q.K=a[1]; q.n=a[2]; q.tandelta=tan(deg2rad(a[3])); }
    const double g = 9.81, gx = g*sin(th);
    basal_state s; s.h=h; s.rho=rho; s.gn=g*cos(th); s.g=g; s.sigma_n=rho*s.gn*h;
    const double dt = 1.0e-6; double u = 0.0, t = 0.0, next = 0.0;
    while(t<=tend+1e-12)
    {
        if(t>=next-1e-12) { printf("%.4f %.10f\n",t,u); next += 0.05; }
        u = implicit_speed(q,s,u + dt*gx,dt);
        t += dt;
    }
    return 0;
}
