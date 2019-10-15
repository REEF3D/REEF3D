/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"increment.h"

#ifndef WAVE_LIB_WCP_H_
#define WAVE_LIB_WCP_H_

using namespace std;

class wave_lib_wcp : public wave_lib_precalc, public wave_lib_parameters, public increment
{
public:
    wave_lib_wcp(lexer*, ghostcell*);
	virtual ~wave_lib_wcp();
    
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
    
    void filename_header(lexer*,ghostcell*);
    void filename(lexer*,ghostcell*,int);
    
    void allocate(lexer*,ghostcell*);
    
    void time_interpol(lexer*);
    void sigma_update(lexer*);
    
    // interpolation
    double space_interpol(lexer*,double,double,double);
    double plane_interpol(lexer*,double,double);
    
    double Upol(lexer*,double,double,double);
    double Vpol(lexer*,double,double,double);
    double Wpol(lexer*,double,double,double);
    double Epol(lexer*,double,double);
    
    int pos_i(double);
    int pos_j(double);
    int pos_k(double);

    // arrays
    double *X;
    double *Y;
    double *Z;
    double **B;
    double *simtime;
    
    double ***U1,***U2,***U;
    double ***V1,***V2,***V;
    double ***W1,***W2,***W;
    double **E1,**E2,**E;
    double **sigz;
    
    
    // variables
    double singamma,cosgamma;
    double val;
    int q1,q2,q1n,q2n;
    double t1,t2,tn,deltaT;
    int Nx,Ny,Nz;
    int num;
    
    int iin;
    float ffn;
	char name[200];
    

    int startup;
    int numiter;
};

#endif
