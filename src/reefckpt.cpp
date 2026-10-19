/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"reefckpt.h"
#include"lexer.h"
#include<mpi.h>
#include<sys/stat.h>
#include<cstdio>
#include<cstring>
#include<algorithm>

using namespace std;

namespace
{
const char MAGIC_PIECE[8] = {'R','3','C','K','P','I','0','1'};
const char MAGIC_MAIN[8]  = {'R','3','C','K','M','N','0','1'};
const int MARGIN = 3;

template<class T> void wr(ofstream &o, const T &v) { o.write((const char*)&v,sizeof(T)); }
template<class T> bool rd(ifstream &i, T &v) { return (bool)i.read((char*)&v,sizeof(T)); }

void wstr(ofstream &o, const string &s)
{
    const int n = (int)s.size();
    wr(o,n);
    o.write(s.data(),n);
}

bool rstr(ifstream &i, string &s)
{
    int n=0;
    if(!rd(i,n) || n<0 || n>4096)
    return false;
    s.resize(n);
    return (bool)i.read(&s[0],n);
}

string piece_name(const string &dir, int r)
{
    char nm[64];
    snprintf(nm,sizeof(nm),"/piece-%06i.dat",r);
    return dir+nm;
}
}

bool reefckpt_open(ofstream &o, const char *name, int I45)
{
    if(I45==1)
    {
        struct stat st;
        if(stat(name,&st)==0)
        {
            o.open(name,ios::app);
            return false;
        }
    }
    o.open(name);
    return true;
}

// ------------------------------------------------------------------------------------------- write

void reefckpt::write_open(lexer *pp, int count)
{
    p = pp;
    char nm[64];
    snprintf(nm,sizeof(nm),"./REEF3D_Checkpoint/%08i",count);
    dirname = nm;
    sc.clear();
    bl.clear();

    if(p->mpirank==0)
    {
        mkdir("./REEF3D_Checkpoint",0777);
        mkdir(dirname.c_str(),0777);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    out.open(piece_name(dirname,p->mpirank).c_str(),ios::binary);
    out.write(MAGIC_PIECE,8);
    const int hd[8] = {p->mpirank,p->origin_i,p->origin_j,p->origin_k,p->knox,p->knoy,p->knoz,MARGIN};
    for(int v : hd)
    wr(out,v);
}

long long reefckpt::count_of(int type) const
{
    const long long s = (long long)p->imax*p->jmax;
    if(type==TS || type==TSI) return s;
    if(type==TF || type==TFI) return s*p->kmaxF;
    if(type==TC || type==TCI) return s*p->kmax;
    return 0;
}

void reefckpt::put(const string &name, int type, const void *data, long long n, size_t bytes)
{
    wstr(out,name);
    wr(out,type);
    wr(out,n);
    out.write((const char*)data,(streamsize)(n*(long long)bytes));
}

void reefckpt::scalar(const string &name, double v) { sc.push_back({name,v}); }
void reefckpt::blob(const string &name, const vector<double> &v) { bl.push_back({name,v}); }
void reefckpt::rank_blob(const string &name, const vector<double> &v) { put(name,TB,v.data(),(long long)v.size(),sizeof(double)); }
void reefckpt::put_S(const string &name, const double *V)  { put(name,TS,V,count_of(TS),sizeof(double)); }
void reefckpt::put_Si(const string &name, const int *V)    { put(name,TSI,V,count_of(TSI),sizeof(int)); }
void reefckpt::put_F(const string &name, const double *V)  { put(name,TF,V,count_of(TF),sizeof(double)); }
void reefckpt::put_Fi(const string &name, const int *V)    { put(name,TFI,V,count_of(TFI),sizeof(int)); }
void reefckpt::put_C(const string &name, const double *V)  { put(name,TC,V,count_of(TC),sizeof(double)); }
void reefckpt::put_Ci(const string &name, const int *V)    { put(name,TCI,V,count_of(TCI),sizeof(int)); }

void reefckpt::write_close(lexer *pp)
{
    p = pp;
    const string endtag = "END";
    wstr(out,endtag);
    out.close();

    const int np = p->mpi_size;
    int mybox[6] = {p->origin_i,p->origin_j,p->origin_k,p->knox,p->knoy,p->knoz};
    vector<int> boxes(6*np);
    MPI_Gather(mybox,6,MPI_INT,boxes.data(),6,MPI_INT,0,MPI_COMM_WORLD);

    if(p->mpirank==0)
    {
        ofstream m((dirname+"/main.dat").c_str(),ios::binary);
        m.write(MAGIC_MAIN,8);
        wr(m,np);
        wr(m,p->gknox); wr(m,p->gknoy); wr(m,p->gknoz);
        for(int v : boxes)
        wr(m,v);
        const int ns = (int)sc.size();
        wr(m,ns);
        for(auto &s : sc)
        {
            wstr(m,s.first);
            wr(m,s.second);
        }
        const int nb = (int)bl.size();
        wr(m,nb);
        for(auto &b : bl)
        {
            wstr(m,b.first);
            const long long n = (long long)b.second.size();
            wr(m,n);
            m.write((const char*)b.second.data(),(streamsize)(n*8));
        }
        m.close();
    }

    // every piece and the main file are complete: the latest checkpoint
    MPI_Barrier(MPI_COMM_WORLD);
    if(p->mpirank==0)
    {
        {
            ofstream l("./REEF3D_Checkpoint/latest.tmp");
            l<<dirname<<endl;
        }
        rename("./REEF3D_Checkpoint/latest.tmp","./REEF3D_Checkpoint/latest");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// -------------------------------------------------------------------------------------------- read

bool reefckpt::read_open(lexer *pp, const string &dir, string &msg)
{
    p = pp;
    dirname = dir;
    if(dirname.empty())
    {
        ifstream l("./REEF3D_Checkpoint/latest");
        if(!l || !(l>>dirname))
        {
            msg = "no checkpoint (REEF3D_Checkpoint/latest missing)";
            return false;
        }
    }

    ifstream m((dirname+"/main.dat").c_str(),ios::binary);
    char mg[8];
    if(!m || !m.read(mg,8) || memcmp(mg,MAGIC_MAIN,8)!=0)
    {
        msg = "checkpoint "+dirname+"/main.dat missing or not a checkpoint";
        return false;
    }
    int np=0;
    rd(m,np); rd(m,GNX); rd(m,GNY); rd(m,GNZ);
    if(GNX!=p->gknox || GNY!=p->gknoy || GNZ!=p->gknoz)
    {
        msg = "checkpoint grid "+to_string(GNX)+" x "+to_string(GNY)+" x "+to_string(GNZ)+" differs from the grid";
        return false;
    }
    box.resize(np);
    for(auto &b : box)
    {
        rd(m,b.oi); rd(m,b.oj); rd(m,b.ok); rd(m,b.ni); rd(m,b.nj); rd(m,b.nk);
        if(b.ok!=0 || b.nk!=GNZ)
        zcut = true;
    }
    if(p->origin_k!=0 || p->knoz!=GNZ)
    zcut = true;
    int ns=0;
    rd(m,ns);
    for(int q=0;q<ns;++q)
    {
        string nm; double v=0.0;
        rstr(m,nm); rd(m,v);
        rsc[nm]=v;
    }
    int nb=0;
    rd(m,nb);
    for(int q=0;q<nb;++q)
    {
        string nm; long long n=0;
        rstr(m,nm); rd(m,n);
        vector<double> v(n);
        m.read((char*)v.data(),(streamsize)(n*8));
        rbl[nm]=v;
    }
    if(!m)
    {
        msg = "checkpoint "+dirname+"/main.dat truncated";
        return false;
    }

    // the same decomposition: every rank finds its own box at its rank
    int ok = (np==p->mpi_size && box[p->mpirank].ok==p->origin_k && box[p->mpirank].oi==p->origin_i && box[p->mpirank].oj==p->origin_j
              && box[p->mpirank].ni==p->knox && box[p->mpirank].nj==p->knoy && box[p->mpirank].nk==p->knoz) ? 1 : 0;
    int all=0;
    MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    same = (all==1);
    mine = same ? p->mpirank : -1;



    // index of every piece file (record offsets)
    pidx.assign(np,{});
    for(int r=0;r<np;++r)
    {
        // on the same decomposition only the own piece is needed (rank blobs: all, read on demand)
        if(same && r!=p->mpirank)
        continue;

        ifstream f(piece_name(dirname,r).c_str(),ios::binary);
        char pm[8];
        if(!f || !f.read(pm,8) || memcmp(pm,MAGIC_PIECE,8)!=0)
        {
            msg = "checkpoint piece "+piece_name(dirname,r)+" missing";
            return false;
        }
        int hd[8];
        for(int &v : hd)
        rd(f,v);
        while(true)
        {
            string nm;
            if(!rstr(f,nm))
            {
                msg = "checkpoint piece "+piece_name(dirname,r)+" truncated";
                return false;
            }
            if(nm=="END")
            break;
            int type=0; long long n=0;
            rd(f,type); rd(f,n);
            const long long off = (long long)f.tellg();
            pidx[r][nm] = {off,n};
            rtype[nm] = type;
            const size_t bytes = (type==TSI || type==TFI || type==TCI) ? sizeof(int) : sizeof(double);
            f.seekg(off + n*(long long)bytes);
        }
    }

    return true;
}

bool reefckpt::has(const string &name) const
{
    return rsc.count(name)>0 || rbl.count(name)>0 || rtype.count(name)>0;
}

double reefckpt::get_scalar(const string &name, double def) const
{
    auto it = rsc.find(name);
    return it==rsc.end() ? def : it->second;
}

vector<double> reefckpt::get_blob(const string &name) const
{
    auto it = rbl.find(name);
    return it==rbl.end() ? vector<double>() : it->second;
}

vector<double> reefckpt::get_rank_blob(const string &name, int r)
{
    vector<double> res;
    if(r<0 || r>=(int)box.size())
    return res;
    ifstream f(piece_name(dirname,r).c_str(),ios::binary);
    char pm[8];
    if(!f || !f.read(pm,8))
    return res;
    int hd[8];
    for(int &v : hd)
    rd(f,v);
    while(true)
    {
        string nm;
        if(!rstr(f,nm) || nm=="END")
        break;
        int type=0; long long n=0;
        rd(f,type); rd(f,n);
        const size_t bytes = (type==TSI || type==TFI || type==TCI) ? sizeof(int) : sizeof(double);
        if(nm==name && type==TB)
        {
            res.resize(n);
            f.read((char*)res.data(),(streamsize)(n*8));
            break;
        }
        f.seekg((long long)f.tellg() + n*(long long)bytes);
    }
    return res;
}

vector<vector<double>> reefckpt::get_rank_blobs(const string &name)
{
    vector<vector<double>> res(box.size());
    for(int r=0;r<(int)box.size();++r)
    res[r] = get_rank_blob(name,r);
    return res;
}

bool reefckpt::get(const string &name, int type, void *V, size_t bytes)
{
    auto t = rtype.find(name);
    int ok = (t!=rtype.end() && t->second==type) ? 1 : 0;
    int all=0;
    MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(!all)
    return false;

    const long long nloc = count_of(type);
    const long long nk = nloc/((long long)p->imax*p->jmax);     // values per column

    // the same decomposition: the own piece as it is
    if(same)
    {
        auto it = pidx[mine].find(name);
        if(it==pidx[mine].end() || it->second.second!=nloc)
        return false;
        ifstream f(piece_name(dirname,mine).c_str(),ios::binary);
        f.seekg(it->second.first);
        f.read((char*)V,(streamsize)(nloc*(long long)bytes));
        return (bool)f;
    }

    // another decomposition that cuts the vertical direction (CFD): cell (I,J,K) from the piece
    // whose interior holds (I,J,K) clamped
    if(zcut && nk>1)
    return get3d(name,V,bytes,nk);

    // another decomposition: column (I,J) from the piece whose interior holds (I,J) clamped
    const int gi0 = p->origin_i, gj0 = p->origin_j;
    vector<char> loaded(box.size(),0);
    vector<vector<char>> buf(box.size());

    auto owner = [&](int I, int J)
    {
        for(int r=0;r<(int)box.size();++r)
        if(I>=box[r].oi && I<box[r].oi+box[r].ni && J>=box[r].oj && J<box[r].oj+box[r].nj)
        return r;
        return -1;
    };

    // pieces in the order of the columns, one owner lookup per distinct clamped column
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
        const int I = ii+gi0, J = jj+gj0;
        const int Ic = min(max(I,0),GNX-1), Jc = min(max(J,0),GNY-1);
        const int r = owner(Ic,Jc);
        if(r<0)
        continue;

        if(!loaded[r])
        {
            auto it = pidx[r].find(name);
            if(it==pidx[r].end())
            return false;
            const long long n = it->second.second;
            const long long pimax = box[r].ni+2*MARGIN, pjmax = box[r].nj+2*MARGIN;
            if(n!=pimax*pjmax*nk)
            return false;
            ifstream f(piece_name(dirname,r).c_str(),ios::binary);
            f.seekg(it->second.first);
            buf[r].resize((size_t)(n*(long long)bytes));
            f.read(buf[r].data(),(streamsize)(n*(long long)bytes));
            loaded[r]=1;
        }

        // local index in the piece r (its own imin = jmin = -MARGIN)
        const long long ri = (I-box[r].oi)+MARGIN, rj = (J-box[r].oj)+MARGIN;
        const long long pjmax = box[r].nj+2*MARGIN;
        if(ri<0 || ri>=box[r].ni+2*MARGIN || rj<0 || rj>=pjmax)
        continue;
        const long long src = (ri*pjmax + rj)*nk;
        const long long dst = ((long long)(ii-p->imin)*p->jmax + (jj-p->jmin))*nk;
        memcpy((char*)V + dst*(long long)bytes, buf[r].data() + src*(long long)bytes, (size_t)(nk*(long long)bytes));
    }
    return true;
}

bool reefckpt::get3d(const string &name, void *V, size_t bytes, long long nk)
{
    const int gi0 = p->origin_i, gj0 = p->origin_j, gk0 = p->origin_k;
    const long long extra = nk - (p->knoz+2*MARGIN);       // 0: cells, 1: nodes
    vector<char> loaded(box.size(),0);
    vector<vector<char>> buf(box.size());

    auto owner = [&](int I, int J, int K)
    {
        for(int r=0;r<(int)box.size();++r)
        if(I>=box[r].oi && I<box[r].oi+box[r].ni && J>=box[r].oj && J<box[r].oj+box[r].nj && K>=box[r].ok && K<box[r].ok+box[r].nk)
        return r;
        return -1;
    };

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    for(long long kk=0; kk<nk; ++kk)
    {
        const int I = ii+gi0, J = jj+gj0, K = int(kk)-MARGIN+gk0;
        const int Ic = min(max(I,0),GNX-1), Jc = min(max(J,0),GNY-1), Kc = min(max(K,0),GNZ-1);
        const int r = owner(Ic,Jc,Kc);
        if(r<0)
        continue;

        const long long pimax = box[r].ni+2*MARGIN, pjmax = box[r].nj+2*MARGIN, pk = box[r].nk+2*MARGIN+extra;
        if(!loaded[r])
        {
            auto it = pidx[r].find(name);
            if(it==pidx[r].end())
            return false;
            const long long n = it->second.second;
            if(n!=pimax*pjmax*pk)
            return false;
            ifstream f(piece_name(dirname,r).c_str(),ios::binary);
            f.seekg(it->second.first);
            buf[r].resize((size_t)(n*(long long)bytes));
            f.read(buf[r].data(),(streamsize)(n*(long long)bytes));
            loaded[r]=1;
        }

        const long long ri = (I-box[r].oi)+MARGIN, rj = (J-box[r].oj)+MARGIN, rk = (K-box[r].ok)+MARGIN;
        if(ri<0 || ri>=pimax || rj<0 || rj>=pjmax || rk<0 || rk>=pk)
        continue;
        const long long src = (ri*pjmax + rj)*pk + rk;
        const long long dst = ((long long)(ii-p->imin)*p->jmax + (jj-p->jmin))*nk + kk;
        memcpy((char*)V + dst*(long long)bytes, buf[r].data() + src*(long long)bytes, bytes);
    }
    return true;
}

bool reefckpt::get_S(const string &name, double *V) { return get(name,TS,V,sizeof(double)); }
bool reefckpt::get_Si(const string &name, int *V)   { return get(name,TSI,V,sizeof(int)); }
bool reefckpt::get_F(const string &name, double *V) { return get(name,TF,V,sizeof(double)); }
bool reefckpt::get_Fi(const string &name, int *V)   { return get(name,TFI,V,sizeof(int)); }
bool reefckpt::get_C(const string &name, double *V) { return get(name,TC,V,sizeof(double)); }
bool reefckpt::get_Ci(const string &name, int *V)   { return get(name,TCI,V,sizeof(int)); }
