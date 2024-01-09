/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include<fstream>
#include"increment.h"

#ifndef WAVE_LIB_HDC_H_
#define WAVE_LIB_HDC_H_

using namespace std;

class wave_lib_hdc : public wave_lib_precalc, public wave_lib_parameters, public increment
{
public:
    wave_lib_hdc(lexer*, ghostcell*);
	virtual ~wave_lib_hdc();
    
    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_fi(lexer*,double,double,double);
    
    virtual void parameters(lexer*,ghostcell*);
    virtual void wave_prestep(lexer*,ghostcell*);
    
private:
    
    // functions
    void read_header(lexer*,ghostcell*);
    void read_result(lexer*,ghostcell*,double**,double***,double***,double***,int);
    void read_result_continuous(lexer*,ghostcell*,double**,double***,double***,double***,int);
    
    void fill_result_continuous(lexer*,ghostcell*);
    
    void filename_header(lexer*,ghostcell*);
    void filename_single(lexer*,ghostcell*,int);
    void filename_continuous(lexer*,ghostcell*);
    
    
    void allocate(lexer*,ghostcell*);
    
    void time_interpol(lexer*);
    
    // interpolation
    double ccpol3D(lexer*,double***,double,double,double);
    double ccpol2D(lexer*,double**,double,double);
    double ccpol2DM(lexer*,double***,double,double);
    double space_interpol(lexer*,double***,double,double,double);
    double plane_interpol(lexer*,double**,double,double);
    
    int pos_i(lexer*,double);
    int pos_j(lexer*,double);
    int pos_k(lexer*,double,int,int);
    
    int ihalf(int,int);

    // arrays
    double *X;
    double *Y;
    double *Zsig;
    double ***Z;
    double **B;
    double *simtime;
    
    double ***U1,***U2,***U;
    double ***V1,***V2,***V;
    double ***W1,***W2,***W;
    double **E1,**E2,**E;
    
    
    
    // variables
    double t_start,t_end;
    int endseries;
    double val;
    int q1,q2,q1n,q2n;
    double t1,t2,tn,deltaT;
    int Nx,Ny,Nz;
    int num;
    
    int iin;
    float ffn;
    double ddn;;
	char name[200];
    ifstream result;
    
    int file_version,file_type;
    int ip1,jp1,kp1;
    int ii,jj,kk;
    int iii,jjj,kkk;
    double xp,yp,zp;
    double wa,wb,wc;
    double v1,v2,v3,v4,v5,v6,v7,v8;
    double x1,x2,x3,x4,y1,y2;
    
    double Xstart,Xend,Ystart,Yend;
    

    int startup;
    int numiter,diter,jdir;
    int file_iter;
};

#endif
