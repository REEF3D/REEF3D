// Standalone test of the active columns of REEFMG (reefmg_core::fine_activity): the loops of the
// solve - residual, products, Krylov updates and reductions, restriction, prolongation, the
// transfers of the device path - run over the columns with an active cell only, and the dry
// columns are cleared once when the set of active columns changes (wetting and drying).
// Architect: Hans Bihs
//
// No MPI library needed: mpi_serial/mpi.h stands in for one rank.
//
// Build (any machine):
//   c++ -O2 -std=c++20 -Impi_serial -I../../src reefmg_dry_columns_test.cpp ../../src/reefmg_core.cpp
//       ../../src/reefmg_device.cpp ../../src/reefmg_device_metal.cpp -o reefmg_dry_columns_test
// Run:
//   ./reefmg_dry_columns_test            checks (exit code 1 on a failure)
//
// Test problem: reefmg_test_problem.h (stretched sigma grid, Dirichlet top, dry columns and a
// bottom step), with the dry block moved: A has its own block, B the block of A shifted, C a
// mostly dry domain (a channel of wet columns).
//
// Checks, for the CPU (fp32 and fp64 storage, solve modes 0, 1, 2) and the reference device
// (N 16 3, N 16 4):
//   1. a solver that solved A, then B (the dry columns change), gives bitwise the solution and
//      the iterations of a new solver that solves B alone; the same for B -> A (columns that were
//      dry become wet) and A -> C -> A
//   2. the solution is finite, the dry columns stay exactly zero in u
//   3. the same rows solved twice (no change of the active columns) give bitwise the same solution

#include "reefmg_core.h"
#include "reefmg_device.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "reefmg_test_problem.h"

// the problem of reefmg_test_problem.h with its dry block replaced: dry columns where dry(i,j)
template<class D>
static problem make_dry(int nx,int ny,int nz,D dry)
{
    problem pr=make_problem(nx,ny,nz,false,false);

    for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
    {
        bool a=!dry(i,j);
        if(nz>2 && i>=2*nx/3 && k<nz/3) a=false;                                   // bottom step
        pr.act[pr.idx(i,j,k)]=a;
    }

    auto active=[&](int i,int j,int k)->bool
    {
        if(i<0 || i>=nx || j<0 || j>=ny || k<0 || k>=nz) return false;
        return pr.act[pr.idx(i,j,k)]!=0;
    };

    rng_state=12345u;
    for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    for(int k=0;k<nz;++k)
    {
        const long q=pr.idx(i,j,k);
        pr.N[q]=pr.S[q]=pr.W[q]=pr.E[q]=pr.T[q]=pr.B[q]=0.0;
        pr.F[q]=0.0;
        if(!pr.act[q]){pr.P[q]=1.0; continue;}

        const double hx=pr.dxn[i+1], hy=pr.dyn[j+1], h=pr.hz[k];
        double n=0,s=0,w=0,e=0,t=0,b=0,dir=0;
        if(active(i+1,j,k)) n=-1.0/(0.5*(hx+pr.dxn[i+2])*hx);
        if(active(i-1,j,k)) s=-1.0/(0.5*(hx+pr.dxn[i  ])*hx);
        if(active(i,j+1,k)) w=-1.0/(0.5*(hy+pr.dyn[j+2])*hy);
        if(active(i,j-1,k)) e=-1.0/(0.5*(hy+pr.dyn[j  ])*hy);
        if(k<nz-1){ if(active(i,j,k+1)) t=-1.0/(0.5*(h+pr.hz[k+1])*h); }
        else      dir=2.0/(h*h);
        if(k>0 && active(i,j,k-1)) b=-1.0/(0.5*(h+pr.hz[k-1])*h);
        pr.N[q]=n; pr.S[q]=s; pr.W[q]=w; pr.E[q]=e; pr.T[q]=t; pr.B[q]=b;
        pr.P[q]=-(n+s+w+e+t+b)+dir;
        pr.F[q]=std::sin(0.2*i)*std::cos(0.15*j)+0.3*rnd();
    }
    return pr;
}

struct config
{
    int device=0;        // N 16
    int bits=32;         // N 14
    int mode=1;          // solve mode
    const char *name="";
};

// a solver that is filled and solved one problem after the other, as REEF3D refills it every step
struct solver
{
    test_op op;
    reefmg_core mg;
    config c;

    solver(const problem &pr,const config &cc) : c(cc)
    {
        mg.set_precision(c.bits);
        mg.set_fine_operator(&op);
        mg.set_device(c.device);
        if(!mg.setup(MPI_COMM_SELF,pr.nx,pr.ny,pr.nz,pr.nx,pr.ny,0,pr.dxn.data(),pr.dyn.data(),0,0))
        {
            std::printf("setup failed: %s\n",mg.err());
            std::exit(2);
        }
    }

    std::vector<double> solve(const problem &pr,int &it)
    {
        op.pr=&pr;
        sc_level &L=mg.fine();
        const long Nc=L.size();
        for(long q=0;q<Nc;++q)
        {
            if(c.bits==32)
            {
                L.pf[q]=(float)pr.P[q]; L.nf[q]=(float)pr.N[q]; L.sf[q]=(float)pr.S[q]; L.wf[q]=(float)pr.W[q];
                L.ef[q]=(float)pr.E[q]; L.tf[q]=(float)pr.T[q]; L.bf[q]=(float)pr.B[q];
            }
            else
            {
                L.p[q]=pr.P[q]; L.n[q]=pr.N[q]; L.s[q]=pr.S[q]; L.w[q]=pr.W[q];
                L.e[q]=pr.E[q]; L.t[q]=pr.T[q]; L.b[q]=pr.B[q];
            }
            L.act[q]=pr.act[q];
            L.f[q]=pr.F[q];
            L.u[q]=0.0;
        }
        mg.coarsen();
        double rel=0.0;
        it=mg.solve_auto(1.0e-10,200,rel,1,1,c.mode);
        return L.u;
    }
};

static bool same(const problem &pr,const std::vector<double> &a,const std::vector<double> &b)
{
    for(int i=0;i<pr.nx;++i)
    for(int j=0;j<pr.ny;++j)
    for(int k=0;k<pr.nz;++k)
    {
        const long q=pr.idx(i,j,k);
        if(std::memcmp(&a[q],&b[q],sizeof(double))!=0) return false;
    }
    return true;
}

static bool dry_zero(const problem &pr,const std::vector<double> &u)
{
    for(int i=0;i<pr.nx;++i)
    for(int j=0;j<pr.ny;++j)
    for(int k=0;k<pr.nz;++k)
    {
        const long q=pr.idx(i,j,k);
        if(!std::isfinite(u[q]) || (!pr.act[q] && u[q]!=0.0)) return false;
    }
    return true;
}

// the sequence of problems seq solved by one solver against a new solver for every problem
static void sequence(const config &c,const std::vector<const problem*> &seq,const char *what)
{
    solver S(*seq[0],c);
    for(size_t n=0;n<seq.size();++n)
    {
        int it1=0, it2=0;
        const std::vector<double> u1=S.solve(*seq[n],it1);
        solver F(*seq[n],c);
        const std::vector<double> u2=F.solve(*seq[n],it2);
        if(n>0)
        check(same(*seq[n],u1,u2) && it1==it2,
              std::string(c.name)+": "+what+", problem "+std::to_string(n)+" after the others = a new solver, bitwise ("
              +std::to_string(it1)+" / "+std::to_string(it2)+" iterations)");
        check(dry_zero(*seq[n],u1),std::string(c.name)+": "+what+", problem "+std::to_string(n)+": the solution finite, the dry cells zero");
    }
}

int main()
{
    const int nx=37, ny=29, nz=8;

    const problem A=make_dry(nx,ny,nz,[&](int i,int j){return i>=nx/3 && i<nx/3+6 && j>=ny/4 && j<ny/4+5;});
    const problem B=make_dry(nx,ny,nz,[&](int i,int j){return i>=nx/3+3 && i<nx/3+9 && j>=ny/4+2 && j<ny/4+7;});
    const problem C=make_dry(nx,ny,nz,[&](int i,int j){return !(j>=ny/2-2 && j<ny/2+3);});   // a wet channel

    std::vector<config> cs;
    {
        config c;
        c.name="CPU fp32 BiCGStab";        c.bits=32; c.mode=1; cs.push_back(c);
        c.name="CPU fp32 V-cycles";        c.bits=32; c.mode=0; cs.push_back(c);
        c.name="CPU fp64 BiCGStab";        c.bits=64; c.mode=1; cs.push_back(c);
        c.name="CPU fp64 V-cycles+fallback"; c.bits=64; c.mode=2; cs.push_back(c);
        c.name="device N 16 3";            c.bits=32; c.mode=1; c.device=3; cs.push_back(c);
        c.name="device N 16 4 (GCR)";      c.bits=32; c.mode=1; c.device=4; cs.push_back(c);
    }

    for(const config &c : cs)
    {
        std::printf("\n%s\n",c.name);
        sequence(c,{&A,&B},"A -> B (the dry block moves)");
        sequence(c,{&B,&A},"B -> A");
        sequence(c,{&A,&C,&A},"A -> C -> A (mostly dry and back)");

        // the same rows twice
        solver S(A,c);
        int i1=0,i2=0;
        const std::vector<double> u1=S.solve(A,i1);
        const std::vector<double> u2=S.solve(A,i2);
        check(same(A,u1,u2) && i1==i2,std::string(c.name)+": the same rows twice, bitwise");
    }

    std::printf("\n%s: %d failure(s)\n",nfail==0? "PASSED" : "FAILED",nfail);
    return nfail==0? 0 : 1;
}
