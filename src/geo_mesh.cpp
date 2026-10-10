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

#include"geo_mesh.h"
#include"gridfile_v2.h"
#include"lexer.h"
#include<iostream>
#include<cstdlib>

geo_mesh::geo_mesh() : tri_x(nullptr), tri_y(nullptr), tri_z(nullptr), ntri(0),
                       solidread(0), toporead(0), geodat(0), dxm(0.0),
                       curv_ni(0), curv_nj(0)
{
}

geo_mesh::~geo_mesh()
{
    release();
}

void geo_mesh::allocate(int num)
{
    release();
    
    ntri = num;
    
    const int n = (num>0) ? num : 1;
    
    tri_x = new double*[n];
    tri_y = new double*[n];
    tri_z = new double*[n];
    
    double *bx = new double[3*n]();
    double *by = new double[3*n]();
    double *bz = new double[3*n]();
    
    for(int q=0; q<n; ++q)
    {
    tri_x[q] = bx + 3*q;
    tri_y[q] = by + 3*q;
    tri_z[q] = bz + 3*q;
    }
}

void geo_mesh::release()
{
    if(tri_x!=nullptr)
    {
    delete [] tri_x[0];
    delete [] tri_y[0];
    delete [] tri_z[0];
    
    delete [] tri_x;
    delete [] tri_y;
    delete [] tri_z;
    }
    
    tri_x = tri_y = tri_z = nullptr;
    ntri = 0;
}

int geo_mesh::count(int role) const
{
    int num=0;
    
    for(const geo_object &ob : obj)
    if(ob.role==role)
    ++num;
    
    return num;
}

int geo_mesh::tricount(int role) const
{
    int num=0;
    
    for(const geo_object &ob : obj)
    if(ob.role==role)
    num += ob.te-ob.ts;
    
    return num;
}

void geo_mesh::read(lexer *p, const char *name)
{
    gridv2::reader rd;
    
    if(!rd.open(name,gridv2::magic_geom))
    {
        if(p->mpirank==0)
        {
        cout<<endl;
        cout<<"!!! "<<rd.error<<" !!!"<<endl;
        cout<<"!!! please regenerate the grid with the current DIVEMesh (grid format v2) !!!"<<endl<<endl;
        }
        exit(1);
    }
    
    gridv2::section sc;
    
    // GHDR
    if(rd.find("GHDR",sc))
    {
        vector<int> iv;
        vector<double> dv;
        sc.get_list(iv,dv);
        
        solidread = iv.size()>0 ? iv[0] : 0;
        toporead  = iv.size()>1 ? iv[1] : 0;
        geodat    = iv.size()>2 ? iv[2] : 0;
        dxm       = dv.size()>0 ? dv[0] : 0.0;
        
        rd.check(sc);
    }
    else
    rd.fail(p,"section GHDR missing");
    
    // OBJS
    obj.clear();
    
    if(rd.find("OBJS",sc))
    {
        const int nobj = sc.get_int();
        
        obj.resize(nobj);
        
        for(geo_object &ob : obj)
        {
            ob.role    = sc.get_int();
            ob.keyword = sc.get_int();
            ob.index   = sc.get_int();
            ob.raymode = sc.get_int();
            ob.invert  = sc.get_int();
            
            const long long t0 = sc.get_i64();
            const long long nt = sc.get_i64();
            
            ob.ts = int(t0);
            ob.te = int(t0+nt);
            
            const int npar = sc.get_int();
            ob.param.resize(npar);
            
            for(double &v : ob.param)
            v = sc.get_double();
        }
        
        rd.check(sc);
    }
    else
    rd.fail(p,"section OBJS missing");
    
    // VERT, TIDX: indexed triangle mesh
    vector<double> vert;
    
    if(rd.find("VERT",sc))
    {
        const long long nv = sc.get_i64();
        
        vert.resize(3*size_t(nv>0 ? nv : 0));
        
        for(double &v : vert)
        v = sc.get_double();
        
        rd.check(sc);
    }
    else
    rd.fail(p,"section VERT missing");
    
    if(rd.find("TIDX",sc))
    {
        const long long nt = sc.get_i64();
        sc.pos = 0;
        
        vector<int> tidx;
        
        if(!sc.get_table(tidx,nt,3))
        rd.fail(p,"section TIDX inconsistent");
        
        allocate(int(nt));
        
        const int nv = int(vert.size()/3);
        
        for(int q=0; q<ntri; ++q)
        for(int v=0; v<3; ++v)
        {
            const int id = tidx[3*q+v];
            
            if(id<0 || id>=nv)
            rd.fail(p,"vertex index out of range");
            
            tri_x[q][v] = vert[3*id];
            tri_y[q][v] = vert[3*id+1];
            tri_z[q][v] = vert[3*id+2];
        }
    }
    else
    rd.fail(p,"section TIDX missing");
    
    // CURV: optional, layout version 1
    curv_ni = curv_nj = 0;
    curv_par.clear(); curv_x.clear(); curv_y.clear(); curv_zb.clear();
    
    if(rd.find("CURV",sc))
    {
        const int nint = sc.get_int();
        vector<int> iv(nint>0 ? nint : 0);
        for(int &v : iv)
        v = sc.get_int();
        
        if(nint<3 || iv[0]!=1 || iv[1]<1 || iv[2]<1)
        rd.fail(p,"section CURV: unknown layout");
        
        const int npar = sc.get_int();
        curv_par.resize(npar>0 ? npar : 0);
        for(double &v : curv_par)
        v = sc.get_double();
        
        const size_t nn = size_t(iv[1]+1)*size_t(iv[2]+1);
        
        if(3*nn*sizeof(double)>sc.size)
        rd.fail(p,"section CURV truncated");
        
        curv_x.resize(nn);
        curv_y.resize(nn);
        curv_zb.resize(nn);
        
        if(!sc.get(curv_x.data(),nn*sizeof(double)) || !sc.get(curv_y.data(),nn*sizeof(double)) || !sc.get(curv_zb.data(),nn*sizeof(double)))
        rd.fail(p,"section CURV truncated");
        
        rd.check(sc);
        
        curv_ni = iv[1];
        curv_nj = iv[2];
    }
    
    if(!rd.ok())
    rd.fail(p,"truncated geometry file");
    
    for(const geo_object &ob : obj)
    if(ob.ts<0 || ob.te>ntri || ob.ts>ob.te)
    rd.fail(p,"inconsistent triangle range");
}
