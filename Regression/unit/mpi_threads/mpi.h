// Threads in place of MPI ranks, for the multi-rank REEFMG unit test
// (reefmg_device_ranks_test.cpp): only what reefmg_core uses.
// Architect: Hans Bihs
//
// Every rank is a thread of one process (mpi_threads::run).  One Cartesian
// communicator (MPI_COMM_WORLD; MPI_Comm_dup and MPI_Comm_split_type return
// it), collectives through a shared slot table and a barrier, point-to-point
// through per-rank mailboxes, shared windows as memory of rank 0.
#ifndef REEF3D_MPI_THREADS_H_
#define REEF3D_MPI_THREADS_H_

#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

typedef int MPI_Comm;
typedef int MPI_Datatype;     // element size in bytes; MPI_INT and MPI_FLOAT both 4: see the op codes
typedef int MPI_Op;
typedef int MPI_Request;
typedef int MPI_Info;
typedef long MPI_Aint;
typedef void* MPI_Win;
struct MPI_Status {int MPI_SOURCE, MPI_TAG, MPI_ERROR;};

#define MPI_COMM_NULL      0
#define MPI_COMM_WORLD     1
#define MPI_COMM_SELF      2
#define MPI_PROC_NULL      (-2)
#define MPI_UNDEFINED      (-32766)
#define MPI_CART           1
#define MPI_DOUBLE         108
#define MPI_FLOAT          104
#define MPI_INT            204
#define MPI_LONG           208
#define MPI_SUM            1
#define MPI_MIN            2
#define MPI_MAX            3
#define MPI_STATUS_IGNORE  ((MPI_Status*)0)
#define MPI_STATUSES_IGNORE ((MPI_Status*)0)
#define MPI_INFO_NULL      0
#define MPI_WIN_NULL       ((MPI_Win)0)
#define MPI_COMM_TYPE_SHARED 1

namespace mpi_threads
{
    inline int tsize(MPI_Datatype t){return t%100;}

    struct message {int src, tag; std::vector<char> data;};
    struct recv     {void *buf; int bytes; int src, tag;};

    struct world
    {
        int np=1, npx=1, npy=1, perx=0, pery=0;
        std::mutex mu;
        std::condition_variable cv;
        int arrived=0;
        long generation=0;
        std::vector<const void*> slot;
        std::vector<std::vector<message> > mailbox;
    };
    inline world W;
    inline thread_local int me=0;
    inline thread_local std::vector<recv> pending;

    inline void barrier()
    {
        std::unique_lock<std::mutex> lk(W.mu);
        const long gen=W.generation;
        if(++W.arrived==W.np)
        {
            W.arrived=0;
            ++W.generation;
            W.cv.notify_all();
        }
        else
        W.cv.wait(lk,[&]{return W.generation!=gen;});
    }

    // every rank publishes a pointer; f runs while all are visible
    inline void with_all(const void *mine,const std::function<void(const std::vector<const void*>&)> &f)
    {
        W.slot[me]=mine;
        barrier();
        f(W.slot);
        barrier();
    }

    // ranks = threads; rank r at (r / npy, r % npy) as MPI_Cart_create orders them
    inline void run(int npx,int npy,int perx,int pery,const std::function<void(int)> &body)
    {
        W.np=npx*npy; W.npx=npx; W.npy=npy; W.perx=perx; W.pery=pery;
        W.arrived=0; W.slot.assign(W.np,0); W.mailbox.assign(W.np,{});
        std::vector<std::thread> th;
        for(int r=0;r<W.np;++r) th.emplace_back([r,&body]{me=r; pending.clear(); body(r);});
        for(auto &t: th) t.join();
    }

    inline int coord_rank(int cx,int cy){return cx*W.npy+cy;}
}

inline int MPI_Comm_rank(MPI_Comm c,int *r){*r=(c==MPI_COMM_SELF)? 0 : mpi_threads::me; return 0;}
inline int MPI_Comm_size(MPI_Comm c,int *s){*s=(c==MPI_COMM_SELF)? 1 : mpi_threads::W.np; return 0;}
inline int MPI_Comm_dup(MPI_Comm c,MPI_Comm *d){*d=c; return 0;}
inline int MPI_Comm_free(MPI_Comm *c){*c=MPI_COMM_NULL; return 0;}
inline int MPI_Comm_split_type(MPI_Comm c,int,int,MPI_Info,MPI_Comm *out){*out=c; return 0;}
inline int MPI_Finalized(int *f){*f=0; return 0;}
inline int MPI_Abort(MPI_Comm,int code){std::fprintf(stderr,"MPI_Abort %d\n",code); std::exit(2);}

inline int MPI_Topo_test(MPI_Comm,int *t){*t=MPI_CART; return 0;}
inline int MPI_Cartdim_get(MPI_Comm,int *nd){*nd=2; return 0;}
inline int MPI_Cart_get(MPI_Comm,int,int *dims,int *per,int *crd)
{
    using namespace mpi_threads;
    dims[0]=W.npx; dims[1]=W.npy; per[0]=W.perx; per[1]=W.pery;
    crd[0]=me/W.npy; crd[1]=me%W.npy;
    return 0;
}
inline int MPI_Cart_coords(MPI_Comm,int r,int,int *crd)
{
    using namespace mpi_threads;
    crd[0]=r/W.npy; crd[1]=r%W.npy;
    return 0;
}
inline int MPI_Cart_shift(MPI_Comm,int dir,int disp,int *src,int *dst)
{
    using namespace mpi_threads;
    int cx=me/W.npy, cy=me%W.npy;
    const int n=(dir==0)? W.npx : W.npy, per=(dir==0)? W.perx : W.pery, c=(dir==0)? cx : cy;
    auto at=[&](int v)->int
    {
        if(v<0 || v>=n){ if(!per) return MPI_PROC_NULL; v=(v+n)%n; }
        return dir==0? coord_rank(v,cy) : coord_rank(cx,v);
    };
    *src=at(c-disp); *dst=at(c+disp);
    return 0;
}

inline int MPI_Barrier(MPI_Comm c){if(c!=MPI_COMM_SELF) mpi_threads::barrier(); return 0;}

inline int MPI_Bcast(void *buf,int n,MPI_Datatype t,int root,MPI_Comm c)
{
    using namespace mpi_threads;
    if(c==MPI_COMM_SELF) return 0;
    with_all(buf,[&](const std::vector<const void*> &s){ if(me!=root) std::memcpy(buf,s[root],(size_t)n*tsize(t)); });
    return 0;
}

inline int MPI_Allreduce(const void *in,void *out,int n,MPI_Datatype t,MPI_Op op,MPI_Comm c)
{
    using namespace mpi_threads;
    if(c==MPI_COMM_SELF){std::memmove(out,in,(size_t)n*tsize(t)); return 0;}
    std::vector<char> res((size_t)n*tsize(t));
    with_all(in,[&](const std::vector<const void*> &s)
    {
        for(int e=0;e<n;++e)
        {
            if(t==MPI_DOUBLE || t==MPI_FLOAT)
            {
                double v=0.0;
                for(int r=0;r<W.np;++r)
                {
                    const double x=(t==MPI_DOUBLE)? ((const double*)s[r])[e] : ((const float*)s[r])[e];
                    if(r==0) v=x; else if(op==MPI_SUM) v+=x; else if(op==MPI_MIN) v=std::min(v,x); else v=std::max(v,x);
                }
                if(t==MPI_DOUBLE) ((double*)res.data())[e]=v; else ((float*)res.data())[e]=(float)v;
            }
            else
            {
                long v=0;
                for(int r=0;r<W.np;++r)
                {
                    const long x=(t==MPI_INT)? ((const int*)s[r])[e] : ((const long*)s[r])[e];
                    if(r==0) v=x; else if(op==MPI_SUM) v+=x; else if(op==MPI_MIN) v=std::min(v,x); else v=std::max(v,x);
                }
                if(t==MPI_INT) ((int*)res.data())[e]=(int)v; else ((long*)res.data())[e]=v;
            }
        }
    });
    std::memcpy(out,res.data(),res.size());
    return 0;
}

inline int MPI_Allgather(const void *in,int n,MPI_Datatype t,void *out,int,MPI_Datatype,MPI_Comm c)
{
    using namespace mpi_threads;
    const size_t b=(size_t)n*tsize(t);
    if(c==MPI_COMM_SELF){std::memmove(out,in,b); return 0;}
    with_all(in,[&](const std::vector<const void*> &s){ for(int r=0;r<W.np;++r) std::memcpy((char*)out+r*b,s[r],b); });
    return 0;
}

inline int MPI_Allgatherv(const void *in,int,MPI_Datatype t,void *out,const int *cnt,const int *disp,MPI_Datatype,MPI_Comm c)
{
    using namespace mpi_threads;
    const int sz=tsize(t);
    if(c==MPI_COMM_SELF){std::memmove((char*)out+(size_t)disp[0]*sz,in,(size_t)cnt[0]*sz); return 0;}
    with_all(in,[&](const std::vector<const void*> &s)
    { for(int r=0;r<W.np;++r) std::memcpy((char*)out+(size_t)disp[r]*sz,s[r],(size_t)cnt[r]*sz); });
    return 0;
}

inline int MPI_Isend(const void *buf,int n,MPI_Datatype t,int dest,int tag,MPI_Comm c,MPI_Request *req)
{
    using namespace mpi_threads;
    *req=-1;
    if(dest==MPI_PROC_NULL) return 0;
    if(c==MPI_COMM_SELF) dest=me;
    message m; m.src=me; m.tag=tag;
    m.data.assign((const char*)buf,(const char*)buf+(size_t)n*tsize(t));
    std::lock_guard<std::mutex> lk(W.mu);
    W.mailbox[dest].push_back(std::move(m));
    W.cv.notify_all();
    return 0;
}

inline int MPI_Irecv(void *buf,int n,MPI_Datatype t,int src,int tag,MPI_Comm c,MPI_Request *req)
{
    using namespace mpi_threads;
    *req=-1;
    if(src==MPI_PROC_NULL) return 0;
    if(c==MPI_COMM_SELF) src=me;
    pending.push_back({buf,n*tsize(t),src,tag});
    return 0;
}

inline int MPI_Waitall(int,MPI_Request*,MPI_Status*)
{
    using namespace mpi_threads;
    std::unique_lock<std::mutex> lk(W.mu);
    for(recv &R: pending)
    {
        for(;;)
        {
            std::vector<message> &box=W.mailbox[me];
            bool got=false;
            for(size_t m=0;m<box.size();++m)
            if(box[m].src==R.src && box[m].tag==R.tag)
            {
                std::memcpy(R.buf,box[m].data.data(),(size_t)R.bytes);
                box.erase(box.begin()+(long)m);
                got=true;
                break;
            }
            if(got) break;
            W.cv.wait(lk);
        }
    }
    pending.clear();
    return 0;
}

inline int MPI_Sendrecv(const void *sb,int sn,MPI_Datatype st,int dest,int stag,
                        void *rb,int rn,MPI_Datatype rt,int src,int rtag,MPI_Comm c,MPI_Status*)
{
    MPI_Request q[2];
    MPI_Irecv(rb,rn,rt,src,rtag,c,&q[0]);
    MPI_Isend(sb,sn,st,dest,stag,c,&q[1]);
    return MPI_Waitall(2,q,MPI_STATUSES_IGNORE);
}

inline int MPI_Win_allocate_shared(MPI_Aint size,int,MPI_Info,MPI_Comm c,void *baseptr,MPI_Win *win)
{
    using namespace mpi_threads;
    void *p=(size>0)? std::calloc((size_t)size,1) : 0;
    *(void**)baseptr=p;
    if(c==MPI_COMM_SELF){*win=p; return 0;}
    void *root=0;
    with_all(p,[&](const std::vector<const void*> &s){ root=(void*)s[0]; });
    *win=root;
    return 0;
}
inline int MPI_Win_shared_query(MPI_Win win,int,MPI_Aint *size,int *disp,void *baseptr)
{
    *size=0; *disp=1;
    *(void**)baseptr=win;
    return 0;
}
inline int MPI_Win_free(MPI_Win *win)
{
    using namespace mpi_threads;
    barrier();
    if(me==0) std::free(*win);
    *win=MPI_WIN_NULL;
    return 0;
}

#endif
