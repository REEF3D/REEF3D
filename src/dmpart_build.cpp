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

// dmpart::build: the global grid file of DIVEMesh (M 50 1, 2) decomposed for nranks ranks with the
// decomposition of DIVEMesh (dmpart_decomp*.cpp), and the grid file of every rank written in memory
// as print_grid of DIVEMesh writes it (same sections, same order), so read_grid reads it unchanged.

#include"dmpart_decomp.h"
#include"dmpart_api.h"
#include"gridfile_v2.h"
#include<span>
#include<string>

namespace dmpart
{
namespace
{
    // indices 0..num-1 grouped by a key 1..nsub; each group keeps ascending order (as in print_grid)
    struct index_lists
    {
        std::vector<int> start, idx;

        template<class Key>
        void build(int nsub, int num, Key key)
        {
            start.assign(nsub+2,0);

            for(int r=0;r<num;++r)
            {
                const int m=key(r);
                if(m>=1 && m<=nsub)
                ++start[m+1];
            }

            for(int m=1;m<nsub+2;++m)
            start[m]+=start[m-1];

            idx.resize(start[nsub+1]);
            std::vector<int> pos(start.begin(),start.end());

            for(int r=0;r<num;++r)
            {
                const int m=key(r);
                if(m>=1 && m<=nsub)
                idx[pos[m]++]=r;
            }
        }

        std::span<const int> list(int m) const { return {idx.data()+start[m], idx.data()+start[m+1]}; }
        size_t size(int m) const { return size_t(start[m+1]-start[m]); }
    };
}

bool build(const char *fname, int nranks, std::vector<std::vector<char>> &rankbuf, std::string &msg)
{
    gridv2::reader gf;

    if(!gf.open(fname,gridv2::magic_glob))
    {
        msg = gf.error;
        return false;
    }

    gridv2::section sc;

    // ---- GHED
    std::vector<int> iv;
    std::vector<double> dv;
    if(!gf.find("GHED",sc))
    {
        msg = "section GHED missing";
        return false;
    }
    sc.get_list(iv,dv);
    if(iv.size()<33 || dv.size()<16)
    {
        msg = "section GHED inconsistent";
        return false;
    }

    lexer L;
    int ip=0, dp=0;
    L.knox = iv[ip++]; L.knoy = iv[ip++]; L.knoz = iv[ip++];
    const int gxma = iv[ip++], gyma = iv[ip++], gzma = iv[ip++], gmarge = iv[ip++];
    const int C11 = iv[ip++], C12 = iv[ip++], C13 = iv[ip++], C14 = iv[ip++], C15 = iv[ip++], C16 = iv[ip++];
    L.C21 = iv[ip++]; L.C22 = iv[ip++]; L.C23 = iv[ip++];
    const int idir = iv[ip++], jdir = iv[ip++], kdir = iv[ip++];
    const int D10 = iv[ip++], B6 = iv[ip++], geob = iv[ip++];
    L.M11 = iv[ip++]; L.M12 = iv[ip++]; L.M13 = iv[ip++]; L.M20 = iv[ip++];
    L.M30_x = iv[ip++]; L.M30_y = iv[ip++]; L.M30_z = iv[ip++];
    L.M34 = iv[ip++]; L.M35 = iv[ip++];
    ip++;                                   // M 10 of DIVEMesh (not used: the ranks of the run)
    const int nsurf = iv[ip++];

    const double DXM = dv[dp++], DR = dv[dp++], DS = dv[dp++], DT = dv[dp++];
    L.xmin = dv[dp++]; L.ymin = dv[dp++]; L.zmin = dv[dp++];
    L.xmax = dv[dp++]; L.ymax = dv[dp++]; L.zmax = dv[dp++];
    const double gox = dv[dp++], goy = dv[dp++], alpha = dv[dp++];
    L.M21 = dv[dp++]; L.M41 = dv[dp++]; L.M45 = dv[dp++];

    if(gxma!=increment::xma || gyma!=increment::yma || gzma!=increment::zma || gmarge!=increment::marge)
    {
        msg = "margins of the global grid file differ from REEF3D's";
        return false;
    }

    L.M10 = nranks;

    // a fixed number of subdomains per direction (M 30) has to fit the ranks; M 20 5: else chosen
    if(L.M30_x*L.M30_y*L.M30_z!=nranks)
    {
        if(L.M20==3 || L.M20==4)
        {
            msg = "DIVEMesh M 30 "+std::to_string(L.M30_x)+" "+std::to_string(L.M30_y)+" "+std::to_string(L.M30_z)
                + " (M 20 "+std::to_string(L.M20)+") does not give "+std::to_string(nranks)+" ranks";
            return false;
        }
        L.M30_x = L.M30_y = L.M30_z = 0;
    }

    // ---- NODE
    if(!gf.find("NODE",sc))
    {
        msg = "section NODE missing";
        return false;
    }
    L.XNv.resize(L.knox+1+2*increment::marge);
    L.YNv.resize(L.knoy+1+2*increment::marge);
    L.ZNv.resize(L.knoz+1+2*increment::marge);
    for(double &v : L.XNv) v = sc.get_double();
    for(double &v : L.YNv) v = sc.get_double();
    for(double &v : L.ZNv) v = sc.get_double();
    L.XN = L.XNv.data(); L.YN = L.YNv.data(); L.ZN = L.ZNv.data();
    if(!sc.good)
    {
        msg = "section NODE inconsistent";
        return false;
    }

    // ---- MANP, WBOX
    if(L.M20==4)
    {
        if(!gf.find("MANP",sc))
        {
            msg = "section MANP missing (M 20 4)";
            return false;
        }
        L.M31v.resize(L.M30_x+1); L.M32v.resize(L.M30_y+1); L.M33v.resize(L.M30_z+1);
        for(int r=0;r<L.M30_x-1;++r) L.M31v[r] = sc.get_int();
        for(int r=0;r<L.M30_y-1;++r) L.M32v[r] = sc.get_int();
        for(int r=0;r<L.M30_z-1;++r) L.M33v[r] = sc.get_int();
        L.M31 = L.M31v.data(); L.M32 = L.M32v.data(); L.M33 = L.M33v.data();
    }

    for(int c=0;c<7;++c)
    L.M34c[c].assign(std::max(L.M34,1),0.0);
    if(L.M34>0)
    {
        if(!gf.find("WBOX",sc))
        {
            msg = "section WBOX missing (M 34)";
            return false;
        }
        for(int b=0;b<L.M34;++b)
        for(int c=0;c<7;++c)
        L.M34c[c][b] = sc.get_double();
    }
    L.M34_xs = L.M34c[0].data(); L.M34_xe = L.M34c[1].data();
    L.M34_ys = L.M34c[2].data(); L.M34_ye = L.M34c[3].data();
    L.M34_zs = L.M34c[4].data(); L.M34_ze = L.M34c[5].data();
    L.M34_f = L.M34c[6].data();

    // ---- the fields with their margins
    dive A(&L);
    A.i_dir = idir; A.j_dir = jdir; A.k_dir = kdir;

    const size_t n3 = size_t(L.knox+2*increment::xma)*size_t(L.knoy+2*increment::yma)*size_t(L.knoz+2*increment::zma);
    const size_t n2 = size_t(L.knox+2*increment::xma)*size_t(L.knoy+2*increment::yma);

    if(!gf.find("GFLG",sc) || !sc.get_rle(A.flag.feld,n3))
    {
        msg = "section GFLG missing or inconsistent";
        return false;
    }
    if(!gf.find("GSOL",sc) || !sc.get_rle(A.solid.feld,n3))
    {
        msg = "section GSOL missing or inconsistent";
        return false;
    }
    if(!gf.find("GSLF",sc) || !sc.get_rle(A.flagslice.feld,n2))
    {
        msg = "section GSLF missing or inconsistent";
        return false;
    }

    std::vector<int> surf;
    if(!gf.find("SURF",sc) || !sc.get_table(surf,nsurf,5))
    {
        msg = "section SURF missing or inconsistent";
        return false;
    }

    std::vector<double> geobed, dataset;
    if(geob)
    {
        if(!gf.find("GEOB",sc))
        {
            msg = "section GEOB missing";
            return false;
        }
        geobed.resize(size_t(L.knox)*L.knoy);
        for(double &v : geobed) v = sc.get_double();
    }
    if(D10>0)
    {
        if(!gf.find("DATA",sc))
        {
            msg = "section DATA missing";
            return false;
        }
        dataset.resize(size_t(L.knox)*L.knoy);
        for(double &v : dataset) v = sc.get_double();
    }

    // ---- the decomposition of DIVEMesh
    decomp D;
    D.start(&L,&A);

    if(A.mx*A.my*A.mz!=nranks)
    {
        msg = "the decomposition gives "+std::to_string(A.mx*A.my*A.mz)+" subdomains for "+std::to_string(nranks)+" ranks";
        return false;
    }

    // ---- the grid files of the ranks (print_grid of DIVEMesh)
    dive *a = &A;
    lexer *p = &L;
    const int nsub = a->mx*a->my*a->mz;
    const int marge = increment::marge;

    index_lists bsurf, bpara[6], bparaco[6], bparaslice[4], bparacoslice[4];
    bsurf.build(nsub,nsurf,[&](int r){return a->subgrid(surf[5*r],surf[5*r+1],surf[5*r+2]);});

    const int_table *parasf[6] = {&a->para1sf,&a->para2sf,&a->para3sf,&a->para4sf,&a->para5sf,&a->para6sf};
    const int_table *paraco[6] = {&a->para1co,&a->para2co,&a->para3co,&a->para4co,&a->para5co,&a->para6co};
    const int paracount[6] = {a->para1count,a->para2count,a->para3count,a->para4count,a->para5count,a->para6count};
    const int paracocount[6] = {a->paraco1count,a->paraco2count,a->paraco3count,a->paraco4count,a->paraco5count,a->paraco6count};

    for(int d=0;d<6;++d)
    {
        const int_table &sf=*parasf[d], &co=*paraco[d];
        bpara[d].build(nsub,paracount[d],[&](int r){return a->subgrid(sf[r][0],sf[r][1],sf[r][2]);});
        bparaco[d].build(nsub,paracocount[d],[&](int r){return co[r][3];});
    }

    const int_table *paraslicesf[4] = {&a->paraslice1sf,&a->paraslice2sf,&a->paraslice3sf,&a->paraslice4sf};
    const int_table *paracoslicesf[4] = {&a->paracoslice1sf,&a->paracoslice2sf,&a->paracoslice3sf,&a->paracoslice4sf};
    const int paraslicecount[4] = {a->paraslice1count,a->paraslice2count,a->paraslice3count,a->paraslice4count};
    const int paracoslicecount[4] = {a->paracoslice1count,a->paracoslice2count,a->paracoslice3count,a->paracoslice4count};

    for(int d=0;d<4;++d)
    {
        const int_table &sf=*paraslicesf[d], &co=*paracoslicesf[d];
        bparaslice[d].build(nsub,paraslicecount[d],[&](int r){return a->subslice(sf[r][0],sf[r][1]);});
        bparacoslice[d].build(nsub,paracoslicecount[d],[&](int r){return co[r][2];});
    }

    rankbuf.assign(nsub,std::vector<char>());

    for(int aa=1;aa<=a->mx;++aa)
    for(int bb=1;bb<=a->my;++bb)
    for(int cc=1;cc<=a->mz;++cc)
    {
        const int count = ((aa-1)*a->my + (bb-1))*a->mz + cc;

        const int i0 = a->xnode[aa-1], i1 = a->xnode[aa];
        const int j0 = a->ynode[bb-1], j1 = a->ynode[bb];
        const int k0 = a->znode[cc-1], k1 = a->znode[cc];
        const int ni = i1-i0, nj = j1-j0, nk = k1-k0;

        gridv2::writer w;
        size_t pos;

        w.buf.reserve(size_t(ni+nj+nk+6*marge+3)*sizeof(double) + size_t(ni)*size_t(nj)*4*sizeof(double) + 65536);
        w.start(gridv2::magic_grid);

        // HEAD
        pos = w.begin("HEAD");

        const int nint = 62;
        w.put_int(nint);
        w.put_int(p->M10);
        w.put_int(a->subknox[count]);
        w.put_int(a->subknoy[count]);
        w.put_int(a->subknoz[count]);
        w.put_int(a->knox);
        w.put_int(a->knoy);
        w.put_int(a->knoz);
        w.put_int(i0);
        w.put_int(j0);
        w.put_int(k0);
        w.put_int(int(bsurf.size(count)));
        w.put_int(a->para1[count]);
        w.put_int(a->para2[count]);
        w.put_int(a->para3[count]);
        w.put_int(a->para4[count]);
        w.put_int(a->para5[count]);
        w.put_int(a->para6[count]);
        w.put_int(a->paraco1[count]);
        w.put_int(a->paraco2[count]);
        w.put_int(a->paraco3[count]);
        w.put_int(a->paraco4[count]);
        w.put_int(a->paraco5[count]);
        w.put_int(a->paraco6[count]);
        w.put_int(a->paraslice1[count]);
        w.put_int(a->paraslice2[count]);
        w.put_int(a->paraslice3[count]);
        w.put_int(a->paraslice4[count]);
        w.put_int(a->paracoslice1[count]);
        w.put_int(a->paracoslice2[count]);
        w.put_int(a->paracoslice3[count]);
        w.put_int(a->paracoslice4[count]);
        w.put_int(a->nbpara1[count]);
        w.put_int(a->nbpara2[count]);
        w.put_int(a->nbpara3[count]);
        w.put_int(a->nbpara4[count]);
        w.put_int(a->nbpara5[count]);
        w.put_int(a->nbpara6[count]);
        w.put_int(a->mx);
        w.put_int(a->my);
        w.put_int(a->mz);
        w.put_int(C11);
        w.put_int(C12);
        w.put_int(C13);
        w.put_int(C14);
        w.put_int(C15);
        w.put_int(C16);
        w.put_int(p->C21);
        w.put_int(p->C22);
        w.put_int(p->C23);
        w.put_int(a->periodicX[count][0]);
        w.put_int(a->periodicX[count][1]);
        w.put_int(a->periodicX[count][2]);
        w.put_int(a->periodicX[count][3]);
        w.put_int(a->periodicX[count][4]);
        w.put_int(a->periodicX[count][5]);
        w.put_int(a->i_dir);
        w.put_int(a->j_dir);
        w.put_int(a->k_dir);
        w.put_int(D10);
        w.put_int(B6);
        w.put_int(marge);
        w.put_int(count);

        const int ndbl = 22;
        w.put_int(ndbl);
        w.put_double(DXM);
        w.put_double(DR);
        w.put_double(DS);
        w.put_double(DT);
        w.put_double(a->xorig[aa-1]);
        w.put_double(a->yorig[bb-1]);
        w.put_double(a->zorig[cc-1]);
        w.put_double(a->xorig[aa]);
        w.put_double(a->yorig[bb]);
        w.put_double(a->zorig[cc]);
        w.put_double(p->xmin);
        w.put_double(p->ymin);
        w.put_double(p->zmin);
        w.put_double(p->xmax);
        w.put_double(p->ymax);
        w.put_double(p->zmax);
        w.put_double(gox);
        w.put_double(goy);
        w.put_double(alpha);
        w.put_double(0.0);
        w.put_double(0.0);
        w.put_double(0.0);
        w.end(pos);

        // FLAG
        pos = w.begin("FLAG");
        w.put_rle(size_t(ni)*size_t(nj)*size_t(nk),[&](size_t r)
        {
            const int ii = int(r/(size_t(nj)*size_t(nk)));
            const int jj = int((r/size_t(nk))%size_t(nj));
            const int kk = int(r%size_t(nk));
            return a->flag(i0+ii,j0+jj,k0+kk);
        });
        w.end(pos);

        // NODE
        pos = w.begin("NODE");
        for(int ii=i0-marge;ii<=i1+marge;++ii)
        w.put_double(p->XN[ii+marge]);
        for(int jj=j0-marge;jj<=j1+marge;++jj)
        w.put_double(p->YN[jj+marge]);
        for(int kk=k0-marge;kk<=k1+marge;++kk)
        w.put_double(p->ZN[kk+marge]);
        w.end(pos);

        // SURF
        pos = w.begin("SURF");
        {
            const auto lst = bsurf.list(count);
            w.put_table((long long)lst.size(),5,[&](long long r, int c){return surf[5*lst[r]+c] - (c==0 ? i0 : (c==1 ? j0 : (c==2 ? k0 : 0)));});
        }
        w.end(pos);

        // PARA
        pos = w.begin("PARA");
        for(int d=0;d<6;++d)
        {
            const int_table &sf=*parasf[d];
            const auto lst = bpara[d].list(count);
            w.put_table((long long)lst.size(),3,[&](long long r, int c){return sf[lst[r]][c] - (c==0 ? i0 : (c==1 ? j0 : k0));});
        }
        w.end(pos);

        // PACO
        pos = w.begin("PACO");
        for(int d=0;d<6;++d)
        {
            const int_table &co=*paraco[d];
            const auto lst = bparaco[d].list(count);
            w.put_table((long long)lst.size(),3,[&](long long r, int c){return co[lst[r]][c] - (c==0 ? i0 : (c==1 ? j0 : k0));});
        }
        w.end(pos);

        // SLFL
        pos = w.begin("SLFL");
        w.put_rle(size_t(ni)*size_t(nj),[&](size_t r)
        {
            return a->flagslice(i0+int(r/size_t(nj)),j0+int(r%size_t(nj)));
        });
        w.end(pos);

        // SLPA
        pos = w.begin("SLPA");
        for(int d=0;d<4;++d)
        {
            const int_table &sf=*paraslicesf[d];
            const auto lst = bparaslice[d].list(count);
            w.put_table((long long)lst.size(),2,[&](long long r, int c){return sf[lst[r]][c] - (c==0 ? i0 : j0);});
        }
        w.end(pos);

        // SLPC
        pos = w.begin("SLPC");
        for(int d=0;d<4;++d)
        {
            const int_table &co=*paracoslicesf[d];
            const auto lst = bparacoslice[d].list(count);
            w.put_table((long long)lst.size(),2,[&](long long r, int c){return co[lst[r]][c] - (c==0 ? i0 : j0);});
        }
        w.end(pos);

        // GEOB
        if(geob)
        {
            pos = w.begin("GEOB");
            for(int ii=i0;ii<i1;++ii)
            for(int jj=j0;jj<j1;++jj)
            w.put_double(geobed[size_t(ii)*L.knoy+jj]);
            w.end(pos);
        }

        // DATA
        if(D10>0)
        {
            pos = w.begin("DATA");
            for(int ii=i0;ii<i1;++ii)
            for(int jj=j0;jj<j1;++jj)
            w.put_double(dataset[size_t(ii)*L.knoy+jj]);
            w.end(pos);
        }

        w.finish();
        rankbuf[count-1] = std::move(w.buf);
    }

    return true;
}
}
