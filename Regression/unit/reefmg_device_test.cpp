// Standalone test of the REEFMG device path (N 16): the V-cycle kernels of
// src/reefmg_device_kernels.h on the reference device (N 16 2, any machine)
// and on the Apple Silicon GPU (N 16 1, macOS with -DREEF3D_USE_METAL),
// against the CPU V-cycle of reefmg_core.
// Architect: Hans Bihs
//
// No MPI library needed: mpi_serial/mpi.h stands in for one rank (and its
// periodic self-exchange).
//
// Build (any machine, reference device):
//   g++ -O2 -std=c++20 -Impi_serial -I../../src reefmg_device_test.cpp
//       ../../src/reefmg_core.cpp ../../src/reefmg_device.cpp -o reefmg_device_test
// Build (macOS, reference device and Metal GPU):
//   clang++ -O2 -std=c++20 -fobjc-arc -DREEF3D_USE_METAL -Impi_serial -I../../src reefmg_device_test.cpp
//       ../../src/reefmg_core.cpp ../../src/reefmg_device.cpp ../../src/reefmg_device_metal.mm
//       -framework Metal -framework Foundation -o reefmg_device_test
// Run:
//   ./reefmg_device_test                 checks (exit code 1 on a failure)
//   ./reefmg_device_test bench NX NY NZ  timing of the solvers on one problem
//
// Test problem: the 7-point Poisson operator of a sigma-like grid, stretched
// in x and z, Neumann at the sides and the bottom, Dirichlet at the top (the
// free surface), a block of dry columns and a bottom step of inactive cells
// (identity rows), right-hand side a smooth field plus noise.  Variants: a
// single layer (nz = 1, the depth-averaged solvers), periodic in x (halo
// exchange with oneself, which makes the device synchronise before every
// half-sweep), fp64 storage, a truncated hierarchy.
//
// Checks:
//   1. one V-cycle from u = 0 on the device = the CPU red-black V-cycle (N 12 2)
//      to single-precision round-off (the CPU sums in double)
//   2. BiCGStab (N 10 1 default mode) to 1e-10: the device needs at most one
//      iteration more than the CPU red-black cycle, and its solution equals
//      the CPU fp64 lexicographic solution to 1e-8
//   3. V-cycles alone (mode 0) and with fallback (mode 2) converge on the device
//   4. Metal (if built): one V-cycle = the reference device to 1e-5, the
//      same BiCGStab iterations (+-1)

#include "reefmg_core.h"
#include "reefmg_device.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "reefmg_test_problem.h"

// ------------------------------------------------------------------ one run
struct config
{
    int device=0;        // N 16
    int ordering=0;      // N 12 2 -> 1
    int bits=32;         // N 14
    int mode=1;          // solve mode
    int maxiter=200;
    double tol=1.0e-10;
    int levels=0;        // >0: set_active_levels
    int dcycles=0;       // N 17
    int dsweeps=0;       // N 19
};

struct result
{
    std::vector<double> u;
    int it=0;
    double relres=0.0, seconds=0.0;
    double pcseconds=0.0;    // in the preconditioner
    long pccalls=0;
    std::string timing, profile;
    bool on_device=false;
    std::string status;
};

static result run(const problem &pr,const config &c,int repeat=1)
{
    result res;
    test_op op;
    op.pr=&pr;

    reefmg_core mg;
    mg.set_precision(c.bits);
    mg.set_ordering(c.ordering);
    mg.set_fine_operator(&op);
    mg.set_device(c.device);
    mg.set_device_cycle(c.dcycles,c.dsweeps);

    if(!mg.setup(MPI_COMM_SELF,pr.nx,pr.ny,pr.nz,pr.nx,pr.ny,0,pr.dxn.data(),pr.dyn.data(),pr.perx? 1 : 0,0))
    {
        std::printf("setup failed: %s\n",mg.err());
        std::exit(2);
    }
    if(c.levels>0) mg.set_active_levels(c.levels);
    res.on_device=mg.on_device();
    res.status=mg.device_status();

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
    }

    mg.coarsen();

    double best=1.0e30;
    for(int r=0;r<repeat;++r)
    {
        for(long q=0;q<Nc;++q){L.f[q]=pr.F[q]; L.u[q]=0.0;}
        mg.reset_timers();
        const auto t0=std::chrono::steady_clock::now();
        res.it=mg.solve_auto(c.tol,c.maxiter,res.relres,1,1,c.mode);
        const auto t1=std::chrono::steady_clock::now();
        const double t=std::chrono::duration<double>(t1-t0).count();
        if(t<best)
        {
            best=t;
            res.pcseconds=mg.precondition_time();
            res.pccalls=mg.precondition_calls();
            res.timing=mg.device_timing();
        }
    }
    res.seconds=best;
    res.profile=mg.device_profile();
    res.u=L.u;
    return res;
}

static double maxdiff(const problem &pr,const std::vector<double> &a,const std::vector<double> &b,double *scale=0)
{
    double d=0.0, m=0.0;
    for(int i=0;i<pr.nx;++i)
    for(int j=0;j<pr.ny;++j)
    for(int k=0;k<pr.nz;++k)
    {
        const long q=pr.idx(i,j,k);
        d=std::max(d,std::fabs(a[q]-b[q]));
        m=std::max(m,std::fabs(a[q]));
    }
    if(scale) *scale=m;
    return m>0.0? d/m : d;
}

// ------------------------------------------------------------------ checks
static void suite(const char *name,const problem &pr,int bits,int levels,bool metal)
{
    std::printf("\n%s: %d x %d x %d%s, fp%d%s\n",name,pr.nx,pr.ny,pr.nz,pr.perx? ", periodic x" : "",bits,
                levels>0? ", truncated hierarchy" : "");

    config ref;  ref.bits=64;                         // CPU fp64 lexicographic, the reference solution
    config zeb;  zeb.bits=bits; zeb.ordering=1;       // CPU red-black, the algorithm the device runs
    config dev;  dev.bits=bits; dev.device=2;         // reference device
    ref.levels=zeb.levels=dev.levels=levels;

    // 1. one V-cycle (mode 0, one iteration: defect correction from u = 0)
    {
        config z1=zeb, d1=dev; z1.mode=d1.mode=0; z1.maxiter=d1.maxiter=1; z1.tol=d1.tol=0.0;
        const result rz=run(pr,z1), rd=run(pr,d1);
        check(rd.on_device,"reference device in use ("+rd.status+")");
        const double e=maxdiff(pr,rz.u,rd.u);
        check(e<2.0e-5,fmt("one V-cycle: device = CPU red-black to %.1e (max relative difference)",e));

        if(metal)
        {
            config m1=d1; m1.device=1;
            const result rm=run(pr,m1);
            check(rm.on_device,"Metal device in use ("+rm.status+")");
            const double em=maxdiff(pr,rd.u,rm.u);
            check(em<1.0e-5,fmt("one V-cycle: Metal = reference device to %.1e",em));
        }
    }

    // 2. BiCGStab to 1e-10
    const result r0=run(pr,ref), rz=run(pr,zeb), rd=run(pr,dev);
    check(r0.relres<=1.0e-10 && rz.relres<=1.0e-10 && rd.relres<=1.0e-10,
          fmt("BiCGStab converged: CPU fp64 %.1e, CPU red-black %.1e, device %.1e",r0.relres,rz.relres,rd.relres));
    check(rd.it<=rz.it+1,fmt("BiCGStab iterations: device %g, CPU red-black %g (CPU lexicographic fp64 %g)",
                             rd.it,rz.it,r0.it));
    const double es=maxdiff(pr,r0.u,rd.u);
    check(es<1.0e-8,fmt("device solution = CPU fp64 solution to %.1e",es));

    if(metal)
    {
        config m=dev; m.device=1;
        const result rm=run(pr,m);
        check(rm.relres<=1.0e-10 && std::abs(rm.it-rd.it)<=1,
              fmt("Metal BiCGStab: %g iterations (reference device %g), relres %.1e",rm.it,rd.it,rm.relres));
        const double em=maxdiff(pr,r0.u,rm.u);
        check(em<1.0e-8,fmt("Metal solution = CPU fp64 solution to %.1e",em));
    }

    // stronger device preconditioner (N 17 2, N 19 2): no more iterations, same solution
    {
        config d2=dev; d2.dcycles=2; d2.dsweeps=2;
        const result rs=run(pr,d2);
        const double e2=maxdiff(pr,r0.u,rs.u);
        check(rs.relres<=1.0e-10 && rs.it<=rd.it && e2<1.0e-8,
              fmt("2 V-cycles with 2 sweeps (N 17 2, N 19 2): %g iterations (one cycle %g), = CPU solution to %.1e",rs.it,rd.it,e2));
    }

    // 3. stationary V-cycles (defect correction) and the fallback mode
    {
        config d0=dev, d2=dev; d0.mode=0; d2.mode=2; d0.tol=d2.tol=1.0e-8; d0.maxiter=d2.maxiter=400;
        const result a=run(pr,d0), b=run(pr,d2);
        check(a.relres<=1.0e-8,fmt("device V-cycles (mode 0): %g cycles to %.1e",a.it,a.relres));
        check(b.relres<=1.0e-8,fmt("device V-cycles with fallback (mode 2): %g iterations to %.1e",b.it,b.relres));
    }
}

static void bench(int nx,int ny,int nz)
{
    const problem pr=make_problem(nx,ny,nz,false,true);
    std::printf("benchmark %d x %d x %d = %.2f M cells, BiCGStab to 1e-8, best of 3\n",nx,ny,nz,1.0e-6*nx*ny*nz);

    struct row {const char *name; config c;};
    std::vector<row> rows;
    {config c; c.tol=1e-8; rows.push_back({"CPU fp32 lexicographic (default)",c});}
    {config c; c.tol=1e-8; c.ordering=1; rows.push_back({"CPU fp32 red-black (N 12 2)",c});}
#ifdef REEF3D_USE_METAL
    {config c; c.tol=1e-8; c.device=1; rows.push_back({"Metal GPU (N 16 1)",c});}
    {config c; c.tol=1e-8; c.device=1; c.dsweeps=2; rows.push_back({"Metal, 2 sweeps (N 19 2)",c});}
    {config c; c.tol=1e-8; c.device=1; c.dcycles=2; rows.push_back({"Metal, 2 V-cycles (N 17 2)",c});}
    {config c; c.tol=1e-8; c.device=1; c.dcycles=2; c.dsweeps=2; rows.push_back({"Metal, 2 V-cycles, 2 sweeps",c});}
    {config c; c.tol=1e-8; c.device=1; c.dcycles=3; c.dsweeps=2; rows.push_back({"Metal, 3 V-cycles, 2 sweeps",c});}
#endif
    if(nx*ny*nz<=4000000)
    {config c; c.tol=1e-8; c.device=2; rows.push_back({"reference device (N 16 2)",c});}

    for(size_t r=0;r<rows.size();++r)
    {
        const result res=run(pr,rows[r].c,3);
        std::printf("  %-34s %3d it  %9.4f s  %8.3f ms/it  %7.3f ms/V-cycle (%2.0f %% of the solve)  relres %.1e%s\n",
                    rows[r].name,res.it,res.seconds,1.0e3*res.seconds/std::max(res.it,1),
                    1.0e3*res.pcseconds/std::max(res.pccalls,1L),100.0*res.pcseconds/std::max(res.seconds,1e-30),res.relres,
                    (rows[r].c.device>0 && !res.on_device)? "  (not on the device!)" : "");
        if(rows[r].c.device>0)
        {
            std::printf("  %-34s %s\n","",res.status.c_str());
            std::printf("  %-34s %s\n","",res.timing.c_str());
            if(!res.profile.empty()) std::printf("%s",res.profile.c_str());
        }
    }
}

int main(int argc,char **argv)
{
    if(argc>=5 && std::string(argv[1])=="bench")
    {
        bench(std::atoi(argv[2]),std::atoi(argv[3]),std::atoi(argv[4]));
        return 0;
    }

#ifdef REEF3D_USE_METAL
    bool metal=true;
    {
        std::string err;
        reefmg_device *d=reefmg_device_create(1,err);
        if(!d){std::printf("Metal not available: %s\n",err.c_str()); metal=false;}
        else std::printf("Metal: %s\n",d->name().c_str());
        delete d;
    }
#else
    const bool metal=false;
    std::printf("built without Metal: reference device only\n");
#endif

    suite("3D sigma grid",make_problem(48,40,12,false,true),32,0,metal);
    suite("3D sigma grid, odd extents",make_problem(37,23,9,false,true),32,0,metal);
    suite("3D sigma grid, fp64 storage",make_problem(48,40,12,false,true),64,0,metal);
    suite("3D sigma grid, two levels only",make_problem(32,32,8,false,true),32,2,metal);
    suite("periodic in x",make_problem(64,24,8,true,false),32,0,metal);
    suite("single layer (depth-averaged)",make_problem(96,80,1,false,true),32,0,metal);
    suite("flume, 1 cell wide",make_problem(200,1,10,false,false),32,0,metal);

    // a device that cannot be created leaves the CPU path in place
    {
        problem pr=make_problem(16,16,4,false,false);
        config c; c.device=7;
        const result r=run(pr,c);
        check(!r.on_device && r.relres<=1.0e-10,"unknown device: falls back to the CPU ("+r.status+")");
    }

    std::printf("\n%s: %d failure(s)\n",nfail? "FAILED" : "PASSED",nfail);
    return nfail? 1 : 0;
}
