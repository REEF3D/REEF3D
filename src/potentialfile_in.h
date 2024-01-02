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

#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class turbulence;
class field;

using namespace std;

#ifndef POTENTIALFILE_IN_H_
#define POTENTIALFILE_IN_H_

class potentialfile_in : public increment
{

public:
	potentialfile_in(lexer*,ghostcell*);
	virtual ~potentialfile_in();
    
    virtual void flowfile_start(lexer*,fdm*,ghostcell*,turbulence*);
    virtual void ff_inflow(lexer*,fdm*,ghostcell*,field&,field&,field&);

    virtual void ff_waterlevel(lexer*,fdm*,ghostcell*,field&);

private:
     void filename(lexer*,fdm*,ghostcell*,char*,int);
    void read0(lexer*,fdm*,ghostcell*,turbulence*);
    void read1(lexer*,fdm*,ghostcell*,turbulence*);
     void header_read(lexer*,ghostcell*);
     double ccipol4(lexer*,double**,double,double,double);
     double lint4(double**,int&,int&,int&,double,double,double);
     int conv(double);
     
     ifstream headerfile;
     ifstream potentialfile;

    char name[400];
    char name0[400];
    char name1[400];
    int startup;
    float ffn;
	int iin;
	double ddn;
	int printcount,entrycount;
    int q, qn, count;
    int q0,q1;
    double t0,t1,tn;
    int q0n,q1n;
    double deltaT;
    double deltax;
    int dk,maxk;
    
    double bedlevel,waterlevel;
    
    int Ni,Nj,Nk;
    double xs,xe,ys,ye,zs,ze;
    
    //data
    int *iter;
    double *xloc,*simtime;
    double *T,*S,*E;
    double **U,**W;
     
    
};

#endif



