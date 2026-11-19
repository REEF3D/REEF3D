// DEBUG (NHF_KCHECK): the k index of the 3D macros of iterators3D.h checked against the vertical range of the arrays
// (kmin .. kmin+kmax-1, faces kmaxF); a reach outside is reported once per source line
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
#endif
