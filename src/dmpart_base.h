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

#ifndef DMPART_BASE_H_
#define DMPART_BASE_H_

// dmpart: the domain decomposition of DIVEMesh inside REEF3D (M 40, read_grid).
//
// The dmpart_decomp*.cpp files are the decomposition of DIVEMesh (decomp*.cpp) with the classes
// below in place of the DIVEMesh lexer and dive: only the members the decomposition uses.  They
// run on rank 0 from the global grid file of DIVEMesh (DIVEMesh_Grid/grid-global.dat, DIVEMesh
// M 50 1, 2) and give the same subdomains, parallel surfaces and neighbours as DIVEMesh, so the
// grid of every rank (dmpart_build.cpp) is the one DIVEMesh would write for that number of ranks.
// Only for the dmpart files: the macros of dmpart_looping.h (LOOP, MAX, ...) are DIVEMesh's.

#include"dmpart_looping.h"
#include<vector>
#include<cstddef>
#include<cstdlib>
#include<memory>
#include<new>
#include<iostream>
#include<fstream>
#include<cmath>

namespace dmpart
{
using namespace std;

struct increment
{
protected:
    increment() = default;
    ~increment() = default;

public:
    static inline int i = 0, j = 0, k = 0;
    static inline int aa = 0, bb = 0, cc = 0;
    static inline int n = 0, q = 0, qn = 0;

    static constexpr int xma = 3, yma = 3, zma = 3;
    static constexpr int marge = 5;
};

struct resize_class
{
    template<class T>
    static void arr1(T *&f, int ni)
    {
        if(ni>0)
        f = new T[ni];
        for(int a=0; a<ni; ++a)
        f[a] = T(0);
    }

    void Iarray(int *&f, int ni) { arr1(f,ni); }
    void Darray(double *&f, int ni) { arr1(f,ni); }

    void Iarray(int **&f, int ni, int nj)
    {
        if(ni>0)
        f = new int*[ni];
        for(int a=0; a<ni; ++a)
        arr1(f[a],nj);
    }

    void Darray(double **&f, int ni, int nj)
    {
        if(ni>0)
        f = new double*[ni];
        for(int a=0; a<ni; ++a)
        arr1(f[a],nj);
    }

    void Iarray(int ***&f, int ni, int nj, int nk)
    {
        if(ni>0)
        f = new int**[ni];
        for(int a=0; a<ni; ++a)
        Iarray(f[a],nj,nk);
    }
};

// the members of the DIVEMesh lexer the decomposition uses
class lexer : public increment, public resize_class
{
public:
    int knox=0, knoy=0, knoz=0;
    double xmin=0, ymin=0, zmin=0, xmax=0, ymax=0, zmax=0;
    vector<double> XNv, YNv, ZNv;       // nodes with marge on either side
    double *XN=nullptr, *YN=nullptr, *ZN=nullptr;

    int M10=1, M11=1, M12=1, M13=0, M20=1;
    int M30=0, M30_x=0, M30_y=0, M30_z=0;
    int *M31=nullptr, *M32=nullptr, *M33=nullptr;
    vector<int> M31v, M32v, M33v;
    double M21=0.0, M41=0.1, M45=5.0;
    int M34=0, M35=0;
    vector<double> M34v;                // 7 per box
    double *M34_xs=nullptr,*M34_xe=nullptr,*M34_ys=nullptr,*M34_ye=nullptr,*M34_zs=nullptr,*M34_ze=nullptr,*M34_f=nullptr;
    vector<double> M34c[7];

    int C21=0, C22=0, C23=0;
};

class intfield : public increment
{
public:
    explicit intfield(lexer *p)
    {
        nj = size_t(p->knoy+2*yma);
        nk = size_t(p->knoz+2*zma);
        feld.assign(size_t(p->knox+2*xma)*nj*nk,0);
    }
    intfield(const intfield&) = delete;
    intfield& operator=(const intfield&) = delete;

    inline int& operator()(int ii, int jj, int kk)
    {
        return feld[((size_t)(ii+xma)*nj + (size_t)(jj+yma))*nk + (size_t)(kk+zma)];
    }

    size_t nj, nk;
    vector<int> feld;
};

class intfield2d : public increment
{
public:
    explicit intfield2d(lexer *p)
    {
        nj = size_t(p->knoy+2*yma);
        feld.assign(size_t(p->knox+2*xma)*nj,0);
    }
    intfield2d(const intfield2d&) = delete;
    intfield2d& operator=(const intfield2d&) = delete;

    inline int& operator()(int ii, int jj)
    {
        return feld[(size_t)(ii+xma)*nj + (size_t)(jj+yma)];
    }

    size_t nj;
    vector<int> feld;
};

class int_table
{
public:
    void allocate(size_t rows, int cols)
    {
        ncols = size_t(cols);
        nrows = rows;
        data.reset(static_cast<int*>(std::calloc(rows*ncols>0 ? rows*ncols : 1, sizeof(int))));
        if(!data)
        throw std::bad_alloc();
    }

    inline int* operator[](size_t qq) { return data.get() + qq*ncols; }
    inline const int* operator[](size_t qq) const { return data.get() + qq*ncols; }

    size_t rows() const { return nrows; }

private:
    struct free_deleter { void operator()(int *ptr) const { std::free(ptr); } };
    std::unique_ptr<int,free_deleter> data;
    size_t ncols = 0, nrows = 0;
};

// the members of the DIVEMesh dive the decomposition uses
class dive : public increment, public resize_class
{
public:
    explicit dive(lexer *p) : flag(p), solid(p), subgrid(p), flagslice(p), subslice(p)
    {
        Darray(xorig,p->M10+10);
        Darray(yorig,p->M10+10);
        Darray(zorig,p->M10+10);
        Iarray(xnode,p->M10+10);
        Iarray(ynode,p->M10+10);
        Iarray(znode,p->M10+10);

        knox=p->knox;
        knoy=p->knoy;
        knoz=p->knoz;

        Iarray(subknox,p->M10+1);
        Iarray(subknoy,p->M10+1);
        Iarray(subknoz,p->M10+1);

        Iarray(mpi_index,p->M10+1);
        Iarray(mpi_edges,(p->M10+1)*6);

        Iarray(para1,p->M10+1); Iarray(para2,p->M10+1); Iarray(para3,p->M10+1);
        Iarray(para4,p->M10+1); Iarray(para5,p->M10+1); Iarray(para6,p->M10+1);
        Iarray(paraco1,p->M10+1); Iarray(paraco2,p->M10+1); Iarray(paraco3,p->M10+1);
        Iarray(paraco4,p->M10+1); Iarray(paraco5,p->M10+1); Iarray(paraco6,p->M10+1);
        Iarray(paraslice1,p->M10+1); Iarray(paraslice2,p->M10+1); Iarray(paraslice3,p->M10+1); Iarray(paraslice4,p->M10+1);
        Iarray(paracoslice1,p->M10+1); Iarray(paracoslice2,p->M10+1); Iarray(paracoslice3,p->M10+1); Iarray(paracoslice4,p->M10+1);
        Iarray(nbpara1,p->M10+1); Iarray(nbpara2,p->M10+1); Iarray(nbpara3,p->M10+1);
        Iarray(nbpara4,p->M10+1); Iarray(nbpara5,p->M10+1); Iarray(nbpara6,p->M10+1);
    }

    intfield flag, solid, subgrid;
    intfield2d flagslice, subslice;

    int knox=0, knoy=0, knoz=0;
    int i_dir=1, j_dir=1, k_dir=1;
    int *subknox=nullptr,*subknoy=nullptr,*subknoz=nullptr;

    int *mpx=nullptr,*mpy=nullptr,*mpz=nullptr;
    int mx=1, my=1, mz=1;

    double *xorig=nullptr,*yorig=nullptr,*zorig=nullptr;
    int *xnode=nullptr,*ynode=nullptr,*znode=nullptr;

    int para1count=0,para2count=0,para3count=0,para4count=0,para5count=0,para6count=0;
    int paraslice1count=0,paraslice2count=0,paraslice3count=0,paraslice4count=0;
    int paracoslice1count=0,paracoslice2count=0,paracoslice3count=0,paracoslice4count=0;
    int paraco1count=0,paraco2count=0,paraco3count=0,paraco4count=0,paraco5count=0,paraco6count=0;
    int_table para1sf,para2sf,para3sf,para4sf,para5sf,para6sf;
    int_table paraslice1sf,paraslice2sf,paraslice3sf,paraslice4sf;
    int_table para1co,para2co,para3co,para4co,para5co,para6co;
    int_table paracoslice1sf,paracoslice2sf,paracoslice3sf,paracoslice4sf;
    int *para1,*para2,*para3,*para4,*para5,*para6;
    int *paraco1,*paraco2,*paraco3,*paraco4,*paraco5,*paraco6;
    int *paraslice1,*paraslice2,*paraslice3,*paraslice4;
    int *paracoslice1,*paracoslice2,*paracoslice3,*paracoslice4;
    int *nbpara1,*nbpara2,*nbpara3,*nbpara4,*nbpara5,*nbpara6;

    int*** sgfield=nullptr;
    int **periodicX=nullptr,*periodicXall=nullptr;
    int *mpi_index=nullptr, *mpi_edges=nullptr, mpi_edgenum=0;
};
}

#endif
