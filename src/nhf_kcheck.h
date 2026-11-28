// DEBUG (NHF_KCHECK): the k index of the 3D macros of iterators3D.h checked against the vertical range of the arrays
// (kmin .. kmin+kmax-1, faces kmaxF); a reach outside is reported once per source line. The j index of the macros of
// iterators2D.h and iterators3D.h against jmin .. jmin+jmax-1 likewise (NHF_JCHECK), and the rows of the slice accessor
// (slice_base.h) once per caller address (NHF_JSLICE; addr2line -f -C -e REEF3D <address-1> names the function).
// Build: make all with -DNHF_KCHECK added to its CXXFLAGS
#ifndef NHF_KCHECK_H_
#define NHF_KCHECK_H_
#include<cstdio>
#include<cstdlib>
#include<set>
#include<string>
#include<mutex>
template<class P>
inline int kchk(const P *p, int kk, int face, const char *file, int line)
{
    const int km = face ? p->kmaxF : p->kmax;
    if(p->kmin>-3 && (kk<p->kmin || kk>=p->kmin+km))
    {
        static std::set<std::string> seen;
        char b[512]; snprintf(b,512,"%s:%d",file,line);
        if(seen.insert(b).second)
        fprintf(stderr,"NHF_KCHECK %s k=%d range [%d,%d)%s\n",b,kk,p->kmin,p->kmin+km,face?" F":"");
    }
    return kk;
}
template<class P>
inline int jchk(const P *p, int jj, const char *file, int line)
{
    if(p->jmin>-3 && (jj<p->jmin || jj>=p->jmin+p->jmax))
    {
        static std::set<std::string> seen;
        char b[512]; snprintf(b,512,"%s:%d",file,line);
        if(seen.insert(b).second)
        fprintf(stderr,"NHF_JCHECK %s j=%d range [%d,%d)\n",b,jj,p->jmin,p->jmin+p->jmax);
    }
    return jj;
}
#include<dlfcn.h>
inline void nhf_jslice(void *ra, int jj, int jmin, int jmax)
{
    static std::set<void*> seen;
    if(!seen.insert(ra).second)
    return;
    Dl_info di;
    if(dladdr(ra,&di) && di.dli_fname)
    fprintf(stderr,"NHF_JSLICE %lx %s j=%d range [%d,%d)\n",(unsigned long)((char*)ra-(char*)di.dli_fbase),di.dli_fname,jj,jmin,jmin+jmax);
}
#endif
