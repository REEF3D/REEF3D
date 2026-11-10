// Standalone check of the FNPF sixth-order central differences (A 311 6: fnpf_cds6 and, with
// wetting and drying, fnpf_cds6_wd): no MPI, no REEF3D solver.
// Build:  g++ -O2 -std=c++17 -I../../src fnpf_cds6_test.cpp -o fnpf_cds6_test
// Run:    ./fnpf_cds6_test
// Checks: fx, fy, sx, sy of both classes
//  - give the slope of a linear field exactly on stretched, jittered grids (the form
//    sum(c f)/sum(c x) is exact for linear fields on any grid, if the stencil of the
//    denominator is the one of the numerator),
//  - are exact for y^m, m = 0..6, on a uniform grid (sixth order),
//  - x and y are the same operator: fy of the transposed field on the same grid is bitwise fx.
// The y-stencil of the denominator had YP[JP3] in place of YP[JP2] (-9 YP[JP3]): on a uniform
// grid sum(c y) = 51 dy instead of 60 dy, every d/dy 60/51 = 1.18 times too large.
// Architect: Hans Bihs

// fnpf_cds6*.cpp need only the grid coordinates of the lexer, the indices of increment and
// the element access of field and slice: small stand-ins, so the test builds without the
// solver (the include guards keep the real headers out; looping.h gives the index macros)
#include<cmath>
#include<cstdio>
#include<string>
#include<vector>
#define LEXER_H_
#define FIELD_H_
#define SLICE_H_
#define INCREMENT_H_
#include"looping.h"

class fdm;
class fdm2D;
class increment
{
public:
    increment() {}
    virtual ~increment() {}
    static int i,j,k,n,h;
    static int innercounter;
    static int marge;
    static fdm *aa;
    static fdm2D *bb;
};
int increment::i=0, increment::j=0, increment::k=0, increment::n=0, increment::h=0;
int increment::innercounter=0;
int increment::marge=5;
fdm *increment::aa=nullptr;
fdm2D *increment::bb=nullptr;

class lexer
{
public:
    double *XP=nullptr, *YP=nullptr, *ZP=nullptr, *ZN=nullptr;
    int *flag7=nullptr;
    int knoz=0, A327=0;
    int imin=0, jmin=0, kmin=0, jmax=0, kmaxF=0;
};

static const int NI=16, NJ=16, NK=4, G=5;   // cells and ghost cells per side

class field
{
public:
    field() : v((NI+2*G)*(NJ+2*G)*(NK+2*G),0.0) {}
    double& operator()(int i, int j, int k) {return v[((i+G)*(NJ+2*G) + j+G)*(NK+2*G) + k+G];}
    std::vector<double> v;
};

class slice
{
public:
    slice() : v((NI+2*G)*(NJ+2*G),0.0) {}
    double& operator()(int i, int j) {return v[(i+G)*(NJ+2*G) + j+G];}
    std::vector<double> v;
};

#include"../../src/fnpf_cds6.cpp"
#include"../../src/fnpf_cds6_wd.cpp"

static int nfail = 0, npass = 0;
static void check(bool ok, const std::string &what)
{
    if(ok) ++npass;
    else
    {
        ++nfail;
        if(nfail<=20)
        std::printf("FAIL: %s\n",what.c_str());
    }
}

static bool close(double a, double b, double tol)
{
    return std::fabs(a-b) <= tol*(1.0 + std::fabs(b));
}

// cell centres from nodes: stretched (ratio r) with a deterministic jitter, or uniform
static std::vector<double> centres(int n, double h, double r, double jitter, int seed)
{
    std::vector<double> xn(n+2*G+1), xc(n+2*G);
    double x=-G*h, dx=h;
    unsigned s=seed;
    for(int q=0; q<(int)xn.size(); ++q)
    {
        xn[q]=x;
        s = 1664525u*s + 1013904223u;
        double u = (s>>8)*(1.0/16777216.0) - 0.5;
        x += dx*(1.0 + jitter*u);
        dx *= r;
    }
    for(int q=0; q<(int)xc.size(); ++q)
    xc[q]=0.5*(xn[q]+xn[q+1]);
    return xc;
}

struct grid
{
    std::vector<double> X,Y,Z;
    lexer p;
    grid(double r, double jitter)
    {
        X = centres(NI,0.10,r,jitter,11);
        Y = centres(NJ,0.07,r,jitter,23);
        Z = centres(NK,0.05,1.0,0.0,0);
        p.XP=X.data(); p.YP=Y.data(); p.ZP=Z.data(); p.ZN=Z.data();
    }
    double x(int i) const {return X[i+G];}
    double y(int j) const {return Y[j+G];}
    double z(int k) const {return Z[k+G];}
};

template<class C>
static void linear(const char *name, double r, double jitter)
{
    grid g(r,jitter);
    C d(&g.p);
    const double a=0.3, b=1.7, c=-2.9, e=0.4;
    field f; slice s;
    for(int i=-G; i<NI+G; ++i)
    for(int j=-G; j<NJ+G; ++j)
    {
        s(i,j) = a + b*g.x(i) + c*g.y(j);
        for(int k=-G; k<NK+G; ++k)
        f(i,j,k) = a + b*g.x(i) + c*g.y(j) + e*g.z(k);
    }
    int bad[4]={0,0,0,0};
    for(int i=0; i<NI; ++i)
    for(int j=0; j<NJ; ++j)
    {
        increment::i=i; increment::j=j; increment::k=1;
        if(!close(d.fx(&g.p,f,0.0,0.0),b,1e-12)) ++bad[0];
        if(!close(d.fy(&g.p,f,0.0,0.0),c,1e-12)) ++bad[1];
        if(!close(d.sx(&g.p,s,0.0),b,1e-12)) ++bad[2];
        if(!close(d.sy(&g.p,s,0.0),c,1e-12)) ++bad[3];
    }
    char buf[200];
    const char *op[4]={"fx","fy","sx","sy"};
    for(int q=0; q<4; ++q)
    {
        std::snprintf(buf,sizeof(buf),"%s %s linear field, r=%.2f jitter=%.2f: %d of %d points off",
                      name,op[q],r,jitter,bad[q],NI*NJ);
        check(bad[q]==0,buf);
    }
}

template<class C>
static void sixth_order(const char *name)
{
    grid g(1.0,0.0);
    C d(&g.p);
    for(int m=0; m<=6; ++m)
    {
        field f; slice s;
        for(int i=-G; i<NI+G; ++i)
        for(int j=-G; j<NJ+G; ++j)
        {
            s(i,j) = std::pow(g.y(j),m);
            for(int k=-G; k<NK+G; ++k)
            f(i,j,k) = s(i,j);
        }
        int bad=0;
        for(int i=0; i<NI; ++i)
        for(int j=0; j<NJ; ++j)
        {
            increment::i=i; increment::j=j; increment::k=1;
            double ex = m==0 ? 0.0 : m*std::pow(g.y(j),m-1);
            if(!close(d.fy(&g.p,f,0.0,0.0),ex,1e-9)) ++bad;
            if(!close(d.sy(&g.p,s,0.0),ex,1e-9)) ++bad;
        }
        char buf[200];
        std::snprintf(buf,sizeof(buf),"%s fy/sy of y^%d on a uniform grid: %d of %d points off",name,m,bad,2*NI*NJ);
        check(bad==0,buf);
    }
}

template<class C>
static void symmetry(const char *name)
{
    grid g(1.06,0.3);
    g.Y = g.X;                       // the same coordinates in x and y
    g.p.YP = g.Y.data();
    C d(&g.p);
    field f, ft; slice s, st;
    for(int i=-G; i<NI+G; ++i)
    for(int j=-G; j<NJ+G; ++j)
    {
        double v = std::sin(1.3*g.x(i)) + 0.7*std::cos(2.1*g.x(j)) + 0.2*g.x(i)*g.x(j);
        s(i,j)=v; st(j,i)=v;
        for(int k=-G; k<NK+G; ++k)
        {
            f(i,j,k)=v+0.1*k; ft(j,i,k)=v+0.1*k;
        }
    }
    int bad=0;
    for(int i=0; i<NI; ++i)
    for(int j=0; j<NJ; ++j)
    {
        increment::i=i; increment::j=j; increment::k=1;
        double dfx = d.fx(&g.p,f,0.0,0.0), dsx = d.sx(&g.p,s,0.0);
        increment::i=j; increment::j=i;
        if(d.fy(&g.p,ft,0.0,0.0)!=dfx) ++bad;
        if(d.sy(&g.p,st,0.0)!=dsx) ++bad;
    }
    char buf[200];
    std::snprintf(buf,sizeof(buf),"%s fy/sy of the transposed field bitwise fx/sx: %d of %d points off",name,bad,2*NI*NJ);
    check(bad==0,buf);
}

template<class C>
static void all(const char *name)
{
    linear<C>(name,1.00,0.00);
    linear<C>(name,1.08,0.00);
    linear<C>(name,0.93,0.40);
    sixth_order<C>(name);
    symmetry<C>(name);
}

int main()
{
    all<fnpf_cds6>("fnpf_cds6");
    all<fnpf_cds6_wd>("fnpf_cds6_wd");

    std::printf("fnpf_cds6_test: %d passed, %d failed\n",npass,nfail);
    return nfail==0 ? 0 : 1;
}
