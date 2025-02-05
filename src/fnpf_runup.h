/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs, Edgar Chavez
--------------------------------------------------------------------*/

#include"fieldint5.h"
#include"field5.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;

#ifndef FNPF_RUNUP_H_
#define FNPF_RUNUP_H_

using namespace std;

class fnpf_runup :  public increment
{

public:
	fnpf_runup(lexer*,fdm_fnpf*,ghostcell*,int);
	virtual ~fnpf_runup();
	virtual void start(lexer*,fdm_fnpf*,ghostcell*);
    virtual void ini(lexer*,fdm_fnpf*,ghostcell*);

private:	
	
    void fnpf_runup_calc(lexer*,fdm_fnpf*,ghostcell*);
	void print_fnpf_runup(lexer*,fdm_fnpf*,ghostcell*);
    void print_ini(lexer*,fdm_fnpf*,ghostcell*);
    double acceleration(lexer*, fdm_fnpf*, ghostcell*);
	double dndt(lexer*, fdm_fnpf*, ghostcell*);
	double dudsig(lexer*, fdm_fnpf*, ghostcell*);
	double dvdsig(lexer*, fdm_fnpf*, ghostcell*);
	double dudxi(lexer*, fdm_fnpf*, ghostcell*);
	double dvdxi(lexer*, fdm_fnpf*, ghostcell*);
    double roundFunc(double, int);
	
    // run-up variabes
    double un,vn,xc,yc,rc,etan;
    int k;
    const int ID;
    
    double R1,R2,R3,R4,R5,R6; //Run-up calculations
    
    // printing
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int fnpf_runupprintcount;
    ofstream fout;

    // parallelisation
    double xstart,ystart,xend,yend;
};

#endif


