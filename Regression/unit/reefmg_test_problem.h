// Test problem of the REEFMG unit tests (reefmg_device_test.cpp, reefmg_device_ranks_test.cpp).
// Architect: Hans Bihs
#ifndef REEFMG_TEST_PROBLEM_H_
#define REEFMG_TEST_PROBLEM_H_

#include "reefmg_core.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

inline int nfail=0;
inline void check(bool ok,const std::string &what)
{
    std::printf("%s  %s\n",ok? "  PASS" : "  FAIL",what.c_str());
    if(!ok) ++nfail;
}

inline std::string fmt(const char *f,double a,double b=0.0,double c=0.0)
{
    char s[512];
    std::snprintf(s,sizeof(s),f,a,b,c);
    return s;
}

// ------------------------------------------------------------------ problem
struct problem
{
    int nx,ny,nz;
    bool perx;
    std::vector<double> dxn,dyn,hz;            // widths, dxn/dyn with one ghost each side
    std::vector<double> P,N,S,W,E,T,B,F;       // level layout, halo included
    std::vector<char> act;

    long idx(int i,int j,int k) const {return ((long)(i+1)*(ny+2)+(j+1))*nz+k;}
    long size() const {return (long)(nx+2)*(ny+2)*nz;}
};

inline unsigned rng_state=12345u;
inline double rnd()
{
    rng_state=rng_state*1664525u+1013904223u;
    return (rng_state>>8)*(1.0/16777216.0)-0.5;
}

inline problem make_problem(int nx,int ny,int nz,bool perx,bool holes)
{
    problem pr;
    pr.nx=nx; pr.ny=ny; pr.nz=nz; pr.perx=perx;
    rng_state=12345u;

    pr.dxn.assign(nx+2,1.0); pr.dyn.assign(ny+2,1.0); pr.hz.assign(nz,1.0);
    for(int i=0;i<nx+2;++i) pr.dxn[i]=1.0+0.5*std::sin(0.13*i)*(perx? 0.0 : 1.0);   // periodic: uniform
    for(int j=0;j<ny+2;++j) pr.dyn[j]=1.0;
    double zs=0.0;
    for(int k=0;k<nz;++k){pr.hz[k]=std::pow(1.15,k); zs+=pr.hz[k];}
    for(int k=0;k<nz;++k) pr.hz[k]*=0.3*nz/zs;                                    // stretched sigma layers

    const long Nc=pr.size();
    pr.P.assign(Nc,0.0); pr.N.assign(Nc,0.0); pr.S.assign(Nc,0.0); pr.W.assign(Nc,0.0);
    pr.E.assign(Nc,0.0); pr.T.assign(Nc,0.0); pr.B.assign(Nc,0.0); pr.F.assign(Nc,0.0);
    pr.act.assign(Nc,0);

    // activity: a block of dry columns and a bottom step
    for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
    {
        bool a=true;
        if(holes)
        {
            if(i>=nx/3 && i<nx/3+nx/8+1 && j>=ny/4 && j<ny/4+ny/6+1) a=false;     // dry columns
            if(nz>2 && i>=2*nx/3 && k<nz/3) a=false;                              // bottom step
        }
        pr.act[pr.idx(i,j,k)]=a;
    }

    auto active=[&](int i,int j,int k)->bool
    {
        if(perx){ if(i<0) i+=nx; if(i>=nx) i-=nx; }
        if(i<0 || i>=nx || j<0 || j>=ny || k<0 || k>=nz) return false;
        return pr.act[pr.idx(i,j,k)]!=0;
    };

    for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
    {
        const long q=pr.idx(i,j,k);
        if(!pr.act[q]){pr.P[q]=1.0; continue;}

        const double hx=pr.dxn[i+1], hy=pr.dyn[j+1], h=pr.hz[k];
        double n=0,s=0,w=0,e=0,t=0,b=0,dir=0;
        if(active(i+1,j,k)) n=-1.0/(0.5*(hx+pr.dxn[i+2])*hx);
        if(active(i-1,j,k)) s=-1.0/(0.5*(hx+pr.dxn[i  ])*hx);
        if(active(i,j+1,k)) w=-1.0/(0.5*(hy+pr.dyn[j+2])*hy);
        if(active(i,j-1,k)) e=-1.0/(0.5*(hy+pr.dyn[j  ])*hy);
        if(nz>1)
        {
            if(k<nz-1){ if(active(i,j,k+1)) t=-1.0/(0.5*(h+pr.hz[k+1])*h); }
            else      dir=2.0/(h*h);                                              // free surface, Dirichlet 0
            if(k>0 && active(i,j,k-1)) b=-1.0/(0.5*(h+pr.hz[k-1])*h);
        }
        else
        dir=0.05;                                                                 // depth-averaged: a reaction term
        pr.N[q]=n; pr.S[q]=s; pr.W[q]=w; pr.E[q]=e; pr.T[q]=t; pr.B[q]=b;
        pr.P[q]=-(n+s+w+e+t+b)+dir;
        pr.F[q]=std::sin(0.2*i)*std::cos(0.15*j)+0.3*rnd();
    }
    return pr;
}

// exact fine operator in double, as reefmg::fine_apply provides it from REEF3D's matrix
struct test_op : public sc_operator
{
    const problem *pr;
    void fine_apply(const sc_level &L,const double *x,double *y) override
    {
        const problem &p=*pr;
        const long sx=(long)(L.ny+2)*L.nz, sy=L.nz;
        for(int i=0;i<L.nx;++i)
        for(int j=0;j<L.ny;++j)
        for(int k=0;k<L.nz;++k)
        {
            const long q=L.idx(i,j,k);
            double v=p.P[q]*x[q]+p.N[q]*x[q+sx]+p.S[q]*x[q-sx]+p.W[q]*x[q+sy]+p.E[q]*x[q-sy];
            if(k<L.nz-1) v+=p.T[q]*x[q+1];
            if(k>0)      v+=p.B[q]*x[q-1];
            y[q]=v;
        }
    }
};

#endif
