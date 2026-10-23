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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"part.h"
#include"lexer.h"
#include<cstdint>
#include<cstring>
#include<iostream>

part::part(lexer *p, ghostcell *pgc)
{
    capacity=1;
    index=1;
    index_empty=capacity;
    //ParcelFactor = p->Q41;
    d50 = p->S20;
    rhosed = p->S22;

    ParcelFactor = (1.0-p->S24)*((p->DXM*p->DXM*p->DXM)/(double(MAX(p->Q24,1))*(1.0/6.0)*pow(p->S20,3.0)*PI));

    if(p->mpirank==0)
        cout<<"ParcelFactor: "<<ParcelFactor<<" DXM: "<<p->DXM<<endl;

    //
    p->Darray(U,capacity);
    p->Darray(V,capacity);
    p->Darray(W,capacity);

    p->Darray(URK1,capacity);
    p->Darray(VRK1,capacity);
    p->Darray(WRK1,capacity);

    p->Darray(X,capacity);
    p->Darray(Y,capacity);
    p->Darray(Z,capacity);

    p->Darray(XRK1,capacity);
    p->Darray(YRK1,capacity);
    p->Darray(ZRK1,capacity);

    p->Darray(Uf,capacity);
    p->Darray(Vf,capacity);
    p->Darray(Wf,capacity);

    p->Darray(D,capacity);
    p->Darray(RO,capacity);

    p->Darray(Test,capacity);
    p->Darray(Hop,capacity);
    p->Darray(Wt,capacity);
    Wt[0] = 1.0;

    p->Iarray(Flag,capacity);
    p->Iarray(Empty,capacity);
    
    // Empty[0..index_empty-1] holds the free slots
    Flag[0] = EMPTY;
    Empty[0] = 0;

    // parallel
    capacity_para=1000;
    capacity_buf=1000;

    p->Iarray(sendid,6,capacity_para);

    p->Darray(send,6,capacity_para);
    p->Darray(recv,6,capacity_para);
}

part::~part()
{
    delete[] U;
    U=nullptr;
    delete[] V;
    V=nullptr;
    delete[] W;
    W=nullptr;

    delete[] URK1;
    URK1=nullptr;
    delete[] VRK1;
    VRK1=nullptr;
    delete[] WRK1;
    WRK1=nullptr;

    delete[] X;
    X=nullptr;
    delete[] Y;
    Y=nullptr;
    delete[] Z;
    Z=nullptr;

    delete[] XRK1;
    XRK1=nullptr;
    delete[] YRK1;
    YRK1=nullptr;
    delete[] ZRK1;
    ZRK1=nullptr;

    delete[] Uf;
    Uf=nullptr;
    delete[] Vf;
    Vf=nullptr;
    delete[] Wf;
    Wf=nullptr;

    delete[] D;
    D=nullptr;
    delete[] RO;
    RO=nullptr;

    delete[] Test;
    Test=nullptr;
    delete[] Hop;
    Hop=nullptr;
    delete[] Wt;
    Wt=nullptr;

    delete[] Empty;
    Empty=nullptr;
    delete[] Flag;
    Flag=nullptr;
}

// ---------------------------------------------------------------------------------------------
// checkpoints (P 48, I 45)
// ---------------------------------------------------------------------------------------------

namespace
{
const int PART_ND = 20;     // the double arrays of a slot, in this order
}

#define PART_ARRAYS U,V,W,URK1,VRK1,WRK1,X,Y,Z,XRK1,YRK1,ZRK1,Uf,Vf,Wf,D,RO,Test,Hop,Wt

void part::ckpt_state(std::vector<double> &v) const
{
    for(double x : {double(capacity),double(index),double(index_empty),double(index_empty0),double(numactive),double(numempty),ParcelFactor})
    v.push_back(x);
    for(const double *a : {PART_ARRAYS})
    v.insert(v.end(),a,a+capacity);
    for(int q=0; q<capacity; ++q)
    v.push_back(double(Flag[q]));
    for(int q=0; q<capacity; ++q)
    v.push_back(double(Empty[q]));
}

bool part::ckpt_restore(lexer *p, const std::vector<double> &v, size_t &k)
{
    if(k+7>v.size())
    return false;
    const int cap = int(v[k]);
    if(cap<1 || k+7+size_t(cap)*(PART_ND+2)>v.size())
    return false;
    index = int(v[k+1]); index_empty = int(v[k+2]); index_empty0 = int(v[k+3]);
    numactive = int(v[k+4]); numempty = int(v[k+5]); ParcelFactor = v[k+6];
    k += 7;

    // the arrays at the capacity of the checkpoint (resize decisions depend on it)
    for(double **a : {&U,&V,&W,&URK1,&VRK1,&WRK1,&X,&Y,&Z,&XRK1,&YRK1,&ZRK1,&Uf,&Vf,&Wf,&D,&RO,&Test,&Hop,&Wt})
    {
        delete [] *a;
        *a = new double[cap];
        std::copy(v.begin()+k,v.begin()+k+cap,*a);
        k += cap;
    }
    delete [] Flag;
    delete [] Empty;
    Flag = new int[cap];
    Empty = new int[cap];
    for(int q=0; q<cap; ++q)
    Flag[q] = int(v[k++]);
    for(int q=0; q<cap; ++q)
    Empty[q] = int(v[k++]);
    capacity = cap;

    // the caches of the exchange start again (they are keyed by the position)
    xdest.clear(); xkx.clear(); xky.clear(); xkz.clear(); xki.clear(); xkj.clear(); xkk.clear();
    return true;
}

void part::ckpt_parcels(const std::vector<double> &v, size_t &k, const std::function<bool(double,double,double)> &keep,
                        std::vector<std::vector<double>> &out)
{
    const int cap = int(v[k]);
    const size_t a0 = k+7;
    auto arr = [&](int a, int q){ return v[a0+size_t(a)*cap+q]; };
    const size_t f0 = a0+size_t(PART_ND)*cap;
    for(int q=0; q<cap; ++q)
    {
        const int flag = int(v[f0+q]);
        if(flag==EMPTY || !keep(arr(6,q),arr(7,q),arr(8,q)))
        continue;
        std::vector<double> s(PART_ND+1);
        for(int a=0; a<PART_ND; ++a)
        s[a] = arr(a,q);
        s[PART_ND] = flag;
        out.push_back(s);
    }
    k = f0+2*size_t(cap);
}

void part::ckpt_set(lexer *p, const std::vector<std::vector<double>> &parcels)
{
    // every slot empty, then the parcels in the first slots
    for(int q=0; q<index; ++q)
    Flag[q] = EMPTY;
    index_empty = 0;
    for(int q=index-1; q>=0; --q)
    Empty[index_empty++] = q;
    resize(p,int(parcels.size()));

    // the free slots from the top of the stack are the lowest ones
    numactive = 0;
    for(const std::vector<double> &s : parcels)
    {
        const int q = Empty[--index_empty];
        double *a[PART_ND] = {PART_ARRAYS};
        for(int c=0; c<PART_ND; ++c)
        a[c][q] = s[c];
        Flag[q] = int(s[PART_ND]);
        if(Flag[q]==ACTIVE)
        ++numactive;
    }
    xdest.clear(); xkx.clear(); xky.clear(); xkz.clear(); xki.clear(); xkj.clear(); xkk.clear();
}
#undef PART_ARRAYS
