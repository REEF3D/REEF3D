// Single-rank stand-in for the part of MPI that reefmg_core uses, so that the
// REEFMG unit test (reefmg_device_test.cpp) builds without an MPI library.
// Architect: Hans Bihs
//
// One rank, no Cartesian topology (reefmg_core does not ask for one on a
// single rank).  Point-to-point messages to oneself - the periodic halo of a
// single rank - go through a mailbox: MPI_Isend copies the data in,
// MPI_Waitall delivers it to the matching MPI_Irecv by tag.  Everything to or
// from MPI_PROC_NULL is a no-op.
#ifndef REEF3D_MPI_SERIAL_H_
#define REEF3D_MPI_SERIAL_H_

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>

typedef int MPI_Comm;
typedef int MPI_Datatype;     // the size of one element in bytes
typedef int MPI_Op;
typedef int MPI_Request;      // index into the pending receives, -1 none
struct MPI_Status {int MPI_SOURCE, MPI_TAG, MPI_ERROR;};

#define MPI_COMM_NULL      0
#define MPI_COMM_WORLD     1
#define MPI_COMM_SELF      2
#define MPI_PROC_NULL      (-2)
#define MPI_UNDEFINED      (-32766)
#define MPI_CART           1
#define MPI_DOUBLE         ((MPI_Datatype)sizeof(double))
#define MPI_FLOAT          ((MPI_Datatype)sizeof(float))
#define MPI_INT            ((MPI_Datatype)sizeof(int))
#define MPI_SUM            1
#define MPI_MIN            2
#define MPI_MAX            3
#define MPI_STATUS_IGNORE  ((MPI_Status*)0)
#define MPI_STATUSES_IGNORE ((MPI_Status*)0)
#define MPI_SUCCESS        0

namespace mpi_serial
{
    struct message {int tag; std::vector<char> data;};
    struct recv     {void *buf; int bytes; int tag;};
    inline std::vector<message> mailbox;
    inline std::vector<recv>    pending;
}

inline int MPI_Comm_rank(MPI_Comm,int *r){*r=0; return 0;}
inline int MPI_Comm_size(MPI_Comm,int *s){*s=1; return 0;}
inline int MPI_Comm_dup(MPI_Comm c,MPI_Comm *d){*d=(c==MPI_COMM_NULL? MPI_COMM_SELF : c); return 0;}
inline int MPI_Comm_free(MPI_Comm *c){*c=MPI_COMM_NULL; return 0;}
inline int MPI_Finalized(int *f){*f=0; return 0;}
inline int MPI_Abort(MPI_Comm,int code){std::fprintf(stderr,"MPI_Abort %d\n",code); std::exit(code==0? 1 : 2);}

inline int MPI_Topo_test(MPI_Comm,int *t){*t=MPI_UNDEFINED; return 0;}
inline int MPI_Cartdim_get(MPI_Comm,int *nd){*nd=0; return 0;}
inline int MPI_Cart_get(MPI_Comm,int nd,int *dims,int *per,int *crd)
{
    for(int d=0;d<nd;++d){dims[d]=1; per[d]=0; crd[d]=0;}
    return 0;
}
inline int MPI_Cart_shift(MPI_Comm,int,int,int *a,int *b){*a=*b=MPI_PROC_NULL; return 0;}
inline int MPI_Cart_coords(MPI_Comm,int,int nd,int *crd){for(int d=0;d<nd;++d) crd[d]=0; return 0;}

inline int MPI_Allreduce(const void *in,void *out,int n,MPI_Datatype t,MPI_Op,MPI_Comm)
{
    std::memmove(out,in,(size_t)n*(size_t)t);
    return 0;
}
inline int MPI_Allgather(const void *in,int n,MPI_Datatype t,void *out,int,MPI_Datatype,MPI_Comm)
{
    std::memmove(out,in,(size_t)n*(size_t)t);
    return 0;
}
inline int MPI_Allgatherv(const void *in,int n,MPI_Datatype t,void *out,const int*,const int *disp,MPI_Datatype,MPI_Comm)
{
    std::memmove((char*)out+(size_t)disp[0]*(size_t)t,in,(size_t)n*(size_t)t);
    return 0;
}

inline int MPI_Isend(const void *buf,int n,MPI_Datatype t,int dest,int tag,MPI_Comm,MPI_Request *req)
{
    *req=-1;
    if(dest==MPI_PROC_NULL) return 0;
    mpi_serial::message m;
    m.tag=tag;
    m.data.assign((const char*)buf,(const char*)buf+(size_t)n*(size_t)t);
    mpi_serial::mailbox.push_back(m);
    return 0;
}
inline int MPI_Irecv(void *buf,int n,MPI_Datatype t,int src,int tag,MPI_Comm,MPI_Request *req)
{
    *req=-1;
    if(src==MPI_PROC_NULL) return 0;
    mpi_serial::recv r={buf,n*t,tag};
    mpi_serial::pending.push_back(r);
    *req=(int)mpi_serial::pending.size()-1;
    return 0;
}
inline int MPI_Waitall(int,MPI_Request*,MPI_Status*)
{
    for(size_t r=0;r<mpi_serial::pending.size();++r)
    {
        mpi_serial::recv &R=mpi_serial::pending[r];
        bool found=false;
        for(size_t m=0;m<mpi_serial::mailbox.size();++m)
        if(mpi_serial::mailbox[m].tag==R.tag)
        {
            std::memcpy(R.buf,mpi_serial::mailbox[m].data.data(),(size_t)R.bytes);
            mpi_serial::mailbox.erase(mpi_serial::mailbox.begin()+(long)m);
            found=true;
            break;
        }
        if(!found){std::fprintf(stderr,"mpi_serial: no message for tag %d\n",R.tag); std::exit(2);}
    }
    mpi_serial::pending.clear();
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

#endif
