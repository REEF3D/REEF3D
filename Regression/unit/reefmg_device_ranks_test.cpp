// The REEFMG device path (N 16) on several ranks: threads in place of MPI
// ranks (mpi_threads/mpi.h), so it runs on any machine without MPI.
// Architect: Hans Bihs
//
// Build:
//   c++ -O2 -std=c++20 -pthread -Impi_threads -I../../src reefmg_device_ranks_test.cpp ../../src/reefmg_core.cpp
//       ../../src/reefmg_device.cpp ../../src/reefmg_device_metal.cpp -o reefmg_device_ranks_test
// (on macOS the Metal GPU is checked too)
// Run:  ./reefmg_device_ranks_test
//
// The test problem of reefmg_test_problem.h is cut into the blocks of a
// Cartesian decomposition (2 x 1, 1 x 3, 2 x 2, 3 x 2 with odd block sizes,
// periodic in x on 2 x 2).  Every rank fills its block and solves with
// BiCGStab to 1e-10; the device solves one hierarchy of the whole domain for
// all of them, written and read through the shared window.  Checks:
//   - the device is used on every rank
//   - one V-cycle on np ranks = the same V-cycle on one rank, bitwise: the
//     device sees the same global problem whatever the decomposition
//   - the converged solution = the CPU fp64 solution on one rank to 1e-8
//   - the number of BiCGStab iterations does not depend on the decomposition
//   - two vcycle(0) called by the host (REEFAMR level 0) on np ranks = on one rank, bitwise

#include "reefmg_test_problem.h"
#include "reefmg_device.h"

#include <algorithm>
#include <cstdarg>
#include <mutex>

static std::string sfmt(const char *f,...)
{
    char b[512];
    va_list ap;
    va_start(ap,f);
    std::vsnprintf(b,sizeof(b),f,ap);
    va_end(ap);
    return b;
}

namespace
{
    struct block_result
    {
        std::vector<double> u;   // global interior array (problem layout), this run
        int it=0;
        double relres=0.0;
        bool on_device=false;
        std::string status;
    };

    // problem layout index of a global cell
    long gidx(const problem &g,int i,int j,int k){return g.idx(i,j,k);}

    // the local problem of the block (ox..ox+bnx-1, oy..oy+bny-1), halo included
    problem cut(const problem &g,int ox,int oy,int bnx,int bny)
    {
        problem p;
        p.nx=bnx; p.ny=bny; p.nz=g.nz; p.perx=g.perx;
        p.dxn.assign(bnx+2,1.0); p.dyn.assign(bny+2,1.0); p.hz=g.hz;
        for(int i=-1;i<=bnx;++i)
        {
            int gi=ox+i;
            if(g.perx){ if(gi<0) gi+=g.nx; if(gi>=g.nx) gi-=g.nx; }
            gi=std::max(-1,std::min(g.nx,gi));
            p.dxn[i+1]=g.dxn[gi+1];
        }
        for(int j=-1;j<=bny;++j)
        {
            const int gj=std::max(-1,std::min(g.ny,oy+j));
            p.dyn[j+1]=g.dyn[gj+1];
        }
        const long N=p.size();
        p.P.assign(N,0.0); p.N.assign(N,0.0); p.S.assign(N,0.0); p.W.assign(N,0.0);
        p.E.assign(N,0.0); p.T.assign(N,0.0); p.B.assign(N,0.0); p.F.assign(N,0.0);
        p.act.assign(N,0);
        for(int i=0;i<bnx;++i)
        for(int j=0;j<bny;++j)
        for(int k=0;k<g.nz;++k)
        {
            const long q=p.idx(i,j,k), Q=gidx(g,ox+i,oy+j,k);
            p.P[q]=g.P[Q]; p.N[q]=g.N[Q]; p.S[q]=g.S[Q]; p.W[q]=g.W[Q];
            p.E[q]=g.E[Q]; p.T[q]=g.T[Q]; p.B[q]=g.B[Q]; p.F[q]=g.F[Q];
            p.act[q]=g.act[Q];
        }
        return p;
    }

    // solve the global problem g on an npx x npy decomposition; blocks as
    // DIVEMesh cuts them: the first (n % np) blocks one cell wider
    block_result solve(const problem &g,int npx,int npy,int device,int mode,int maxiter,double tol)
    {
        std::vector<int> bx(npx),by(npy),ox(npx+1,0),oy(npy+1,0);
        for(int c=0;c<npx;++c){bx[c]=g.nx/npx+(c<g.nx%npx?1:0); ox[c+1]=ox[c]+bx[c];}
        for(int c=0;c<npy;++c){by[c]=g.ny/npy+(c<g.ny%npy?1:0); oy[c+1]=oy[c]+by[c];}

        block_result out;
        out.u.assign(g.size(),0.0);
        std::mutex mu;
        bool all_dev=true;

        mpi_threads::run(npx,npy,g.perx? 1 : 0,0,[&](int rank)
        {
            const int cx=rank/npy, cy=rank%npy;
            const problem p=cut(g,ox[cx],oy[cy],bx[cx],by[cy]);
            test_op op; op.pr=&p;

            reefmg_core mg;
            mg.set_precision(32);
            mg.set_fine_operator(&op);
            mg.set_device(device);
            if(!mg.setup(MPI_COMM_WORLD,p.nx,p.ny,p.nz,g.nx,g.ny,0,p.dxn.data(),p.dyn.data(),g.perx? 1 : 0,0))
            {
                std::printf("setup failed: %s\n",mg.err());
                std::exit(2);
            }

            sc_level &L=mg.fine();
            for(long q=0;q<L.size();++q)
            {
                L.pf[q]=(float)p.P[q]; L.nf[q]=(float)p.N[q]; L.sf[q]=(float)p.S[q]; L.wf[q]=(float)p.W[q];
                L.ef[q]=(float)p.E[q]; L.tf[q]=(float)p.T[q]; L.bf[q]=(float)p.B[q];
                L.act[q]=p.act[q]; L.f[q]=p.F[q]; L.u[q]=0.0;
            }
            mg.coarsen();
            double rr=0.0;
            int it=0;
            if(maxiter<0)      // -n: n host-driven vcycle(0) from u = 0 (REEFAMR level 0, reefmg2D)
            for(int n=0;n<-maxiter;++n) mg.vcycle(0,1,1);
            else
            it=mg.solve_auto(tol,maxiter,rr,1,1,mode);

            std::lock_guard<std::mutex> lk(mu);
            for(int i=0;i<p.nx;++i)
            for(int j=0;j<p.ny;++j)
            for(int k=0;k<p.nz;++k)
            out.u[gidx(g,ox[cx]+i,oy[cy]+j,k)]=L.u[L.idx(i,j,k)];
            out.it=it; out.relres=rr;
            all_dev=all_dev && mg.on_device();
            if(rank==0) out.status=mg.device_status();
        });
        out.on_device=all_dev;
        return out;
    }

    double maxrel(const problem &g,const std::vector<double> &a,const std::vector<double> &b)
    {
        double d=0.0,m=0.0;
        for(int i=0;i<g.nx;++i) for(int j=0;j<g.ny;++j) for(int k=0;k<g.nz;++k)
        {
            const long q=g.idx(i,j,k);
            d=std::max(d,std::fabs(a[q]-b[q])); m=std::max(m,std::fabs(a[q]));
        }
        return m>0.0? d/m : d;
    }
}

int main()
{
    std::vector<int> devices={2};
#if REEF3D_METAL
    {
        std::string err;
        reefmg_device *d=reefmg_device_create(1,err);
        if(d){devices.push_back(1); std::printf("Metal: %s\n",d->name().c_str());}
        else std::printf("Metal not available: %s\n",err.c_str());
        delete d;
    }
#endif

    struct dec {int npx,npy;};
    const dec decs[]={{2,1},{1,3},{2,2},{3,2}};

    for(int per=0;per<2;++per)
    {
        const problem g=per? make_problem(40,26,8,true,false) : make_problem(45,31,10,false,true);
        std::printf("\n%s: %d x %d x %d\n",per? "periodic in x" : "sigma grid with dry columns and a step",g.nx,g.ny,g.nz);

        const block_result ref=solve(g,1,1,0,1,200,1.0e-10);            // CPU on one rank (fp32 storage)

        for(int device: devices)
        {
            const block_result one =solve(g,1,1,device,0,1,0.0);       // one V-cycle, one rank
            const block_result host1=solve(g,1,1,device,0,-2,0.0);     // two host-driven vcycle(0)
            const block_result one_bicg=solve(g,1,1,device,1,200,1.0e-10);

            for(const dec &d: decs)
            {
                if(per && d.npx==1) continue;
                char nm[64];
                std::snprintf(nm,sizeof(nm),"%s %d x %d",device==1? "Metal" : "reference",d.npx,d.npy);

                const block_result v=solve(g,d.npx,d.npy,device,0,1,0.0);
                check(v.on_device,std::string(nm)+": device on every rank ("+v.status+")");
                const double ev=maxrel(g,one.u,v.u);
                check(ev==0.0,sfmt("%s: one V-cycle = the one-rank V-cycle, largest difference %.1e",nm,ev));

                const block_result h=solve(g,d.npx,d.npy,device,0,-2,0.0);
                const double eh=maxrel(g,host1.u,h.u);
                check(eh==0.0,sfmt("%s: two host-driven vcycle(0) = the same on one rank, largest difference %.1e",nm,eh));

                const block_result s=solve(g,d.npx,d.npy,device,1,200,1.0e-10);
                const double es=maxrel(g,ref.u,s.u);
                check(s.relres<=1.0e-10 && es<1.0e-8,
                      sfmt("%s: BiCGStab %d iterations, relres %.1e, = CPU solution to %.1e",nm,s.it,s.relres,es));
                check(std::abs(s.it-one_bicg.it)<=1,sfmt("%s: iterations %d, on one rank %d",nm,s.it,one_bicg.it));
            }
        }
    }

    std::printf("\n%s: %d failure(s)\n",nfail? "FAILED" : "PASSED",nfail);
    return nfail? 1 : 0;
}
