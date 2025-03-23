/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include <cstdint>
#include <cstring>
#include<iostream>

part::part(lexer *p, ghostcell *pgc)
{	
    capacity=1;
    index=1;
    index_empty=capacity;
    ParcelFactor = p->Q41;
    d50 = p->S20;
    rhosed = p->S22;
    
    
    ParcelFactor = p->S24*((p->DXM*p->DXM*p->DXM)/(p->Q24*(1.0/6.0)*pow(p->S20,3.0)*PI));
    
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
  
    p->Iarray(Flag,capacity);
    p->Iarray(Empty,capacity);
    
    // parallel
    capacity_para=1000;
    
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
    
    delete[] Empty;
    Empty=nullptr;
    delete[] Flag;
    Flag=nullptr;
}


