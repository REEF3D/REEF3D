// Standalone test of REEFMG with an island of rows without reaction: the coarse diagonal of an
// agglomerate whose fine rows have row sum 0 and are coupled only to each other.
// Architect: Hans Bihs
//
// This is what the wetting and drying iteration of the semi-implicit eta system (si_eta2D, Casulli
// 2009) can give: dry cells (no storage) whose open faces lead only to each other. Before the guard in
// reefmg_core::coarsen_t (and the device kernel) the coarse diagonal of such an agglomerate was 0, the
// line solver divided by it and the V-cycle returned NaN (the storm run of the German Bight case
// stopped after 7.1 h with "HYPRE solver broke down").
//
// No MPI library needed: mpi_serial/mpi.h stands in for one rank.
//
// Build:
//   c++ -O2 -std=c++20 -Impi_serial -I../../src reefmg_island_test.cpp ../../src/reefmg_core.cpp
//       ../../src/reefmg_device.cpp ../../src/reefmg_device_metal.cpp -o reefmg_island_test
// Run:
//   ./reefmg_island_test        checks (exit code 1 on a failure)
//
// Test problem: the depth-averaged screened Poisson operator of reefmg_test_problem.h (nz = 1, reaction
// 0.05) on 16 x 16 cells; a 2 x 2 island on the cells of one coarse agglomerate (i 4-5, j 6-7) and a 1 x 3
// island across two agglomerates (i 11, j 9-11): no reaction, the faces towards the rest closed, the
// right-hand side of each island with zero mean (a consistent singular block).
//
// Checks, CPU fp64 and fp32 storage and the reference device (N 16 2):
//   1. one V-cycle from u = 0 with the islands as they are (no reaction: a singular block) is finite
//      everywhere (before the guard: NaN in fp32)
//   2. the islands with the storage of wet cells (reaction 1, what si_eta2D::wd_islands gives them):
//      BiCGStab with the V-cycle (N 10 1 default mode) converges to 1e-10, and outside the islands its
//      solution equals the solution of the problem with the islands as identity rows to 1e-8
// A singular island itself is not solved: a Krylov method with a preconditioner that does not keep its
// null space diverges on it - si_eta2D keeps such islands out of the system (wd_islands).

#include "reefmg_core.h"
#include "reefmg_device.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "reefmg_test_problem.h"

static bool in_island(int i,int j)
{
    return (i>=4 && i<=5 && j>=6 && j<=7) || (i==11 && j>=9 && j<=11);
}

static int island_id(int i,int j)
{
    if(i>=4 && i<=5 && j>=6 && j<=7) return 1;
    if(i==11 && j>=9 && j<=11) return 2;
    return 0;
}

// islands = true: the islands as described, with the reaction react_isl; false: the island cells as
// identity rows (act 0)
static problem make_island_problem(bool islands,double react_isl=0.0)
{
    const int nx=16, ny=16;
    problem pr=make_problem(nx,ny,1,false,false);

    for(int i=0;i<nx;++i)
    for(int j=0;j<ny;++j)
    {
        const long q=pr.idx(i,j,0);
        const int id=island_id(i,j);
        const double hx=pr.dxn[i+1], hy=pr.dyn[j+1];
        auto cpl=[&](int i2,int j2,double v) -> double
        {
            if(i2<0 || j2<0 || i2>=nx || j2>=ny) return 0.0;
            if(island_id(i2,j2)!=id) return 0.0;          // faces between an island and the rest closed
            return v;
        };
        double n=cpl(i+1,j,-1.0/(0.5*(hx+pr.dxn[i+2])*hx));
        double s=cpl(i-1,j,-1.0/(0.5*(hx+pr.dxn[i  ])*hx));
        double w=cpl(i,j+1,-1.0/(0.5*(hy+pr.dyn[j+2])*hy));
        double e=cpl(i,j-1,-1.0/(0.5*(hy+pr.dyn[j  ])*hy));
        const double react = id ? react_isl : 0.05;

        if(id && !islands)
        {
            pr.act[q]=0;
            pr.P[q]=1.0; pr.N[q]=pr.S[q]=pr.W[q]=pr.E[q]=0.0; pr.F[q]=0.0;
            continue;
        }
        pr.N[q]=n; pr.S[q]=s; pr.W[q]=w; pr.E[q]=e;
        pr.P[q]=-(n+s+w+e)+react;
    }

    // zero-mean right-hand side on each island (consistent)
    for(int id=1;id<=2;++id)
    {
        double sum=0.0; int cnt=0;
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) if(island_id(i,j)==id){sum+=pr.F[pr.idx(i,j,0)]; ++cnt;}
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j) if(island_id(i,j)==id) pr.F[pr.idx(i,j,0)]-=sum/cnt;
    }
    return pr;
}

struct outcome
{
    std::vector<double> u;
    bool finite=true;
    int it=0;
    double relres=1.0;
};

static outcome run(const problem &pr,int bits,int device,bool one_cycle)
{
    outcome o;
    test_op op;
    op.pr=&pr;

    reefmg_core mg;
    mg.set_precision(bits);
    mg.set_fine_operator(&op);
    mg.set_device(device);
    if(!mg.setup(MPI_COMM_SELF,pr.nx,pr.ny,pr.nz,pr.nx,pr.ny,0,pr.dxn.data(),pr.dyn.data(),0,0))
    {
        std::printf("setup failed: %s\n",mg.err());
        std::exit(2);
    }

    sc_level &L=mg.fine();
    const long Nc=L.size();
    for(long q=0;q<Nc;++q)
    {
        if(bits==32)
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
    }

    mg.coarsen();
    for(long q=0;q<Nc;++q){L.f[q]=pr.F[q]; L.u[q]=0.0;}

    if(one_cycle)
    mg.vcycle(0,1,1);
    else
    o.it=mg.solve_auto(1.0e-10,200,o.relres,1,1,1);

    o.u=L.u;
    for(int i=0;i<pr.nx;++i)
    for(int j=0;j<pr.ny;++j)
    if(!std::isfinite(o.u[pr.idx(i,j,0)])) o.finite=false;
    return o;
}

int main()
{
    const problem pi=make_island_problem(true), pw=make_island_problem(true,1.0), pr=make_island_problem(false);

    struct variant {int bits, device; const char *name;};
    for(const variant &v : {variant{64,0,"CPU fp64"}, variant{32,0,"CPU fp32"}, variant{32,2,"reference device (N 16 2)"}})
    {
        std::printf("\n%s\n",v.name);

        const outcome c=run(pi,v.bits,v.device,true);
        check(c.finite,"one V-cycle with the islands is finite");

        const outcome s=run(pw,v.bits,v.device,false);
        check(s.finite && s.relres<=1.0e-10,fmt("BiCGStab with the islands kept wet converges: %.0f iterations, relative residual %.1e",s.it,s.relres));

        const outcome r=run(pr,v.bits,v.device,false);
        double d=0.0, m=0.0;
        for(int i=0;i<pi.nx;++i)
        for(int j=0;j<pi.ny;++j)
        {
            if(in_island(i,j)) continue;
            const long q=pi.idx(i,j,0);
            const double dq=std::fabs(s.u[q]-r.u[q]);
            d=std::isfinite(dq)? std::max(d,dq) : INFINITY;
            m=std::max(m,std::fabs(r.u[q]));
        }
        check(d<=1.0e-8*m,fmt("outside the islands = the problem without them to %.1e (relative)",m>0.0? d/m : d));
    }

    std::printf("\n%s: %d failure(s)\n",nfail? "FAILED" : "passed",nfail);
    return nfail? 1 : 0;
}
