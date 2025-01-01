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
Author: Arun Kamath, Tobias Martin
--------------------------------------------------------------------*/

#include"increment.h"
#include"fieldint5.h"
#include"field5.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;

#ifndef FNPF_PRINT_KINEMATICS_H_
#define FNPF_PRINT_KINEMATICS_H_

using namespace std;

class fnpf_print_kinematics :  public increment
{

public:
	fnpf_print_kinematics(lexer*,fdm_fnpf*,ghostcell*,int);
	virtual ~fnpf_print_kinematics();
	virtual void start(lexer*,fdm_fnpf*,ghostcell*);
    virtual void ini(lexer*,fdm_fnpf*,ghostcell*);

private:	
    void kinematics_calc(lexer*,fdm_fnpf*,ghostcell*);
	void print_kinematics(lexer*,fdm_fnpf*,ghostcell*);
    void print_ini(lexer*,fdm_fnpf*,ghostcell*);
	double dndt(lexer*, fdm_fnpf*, ghostcell*);
	double dudsig(lexer*, fdm_fnpf*, ghostcell*);
	double dvdsig(lexer*, fdm_fnpf*, ghostcell*);
	double dudxi(lexer*, fdm_fnpf*, ghostcell*);
	double dvdxi(lexer*, fdm_fnpf*, ghostcell*);
	
    // force ale variabes
    double Fx1,Fy1,Fx,Fy,xc,yc,rc,cd,cm,etan,dtn,eta2n,ax1,ay1,ax2,ay2,ax3,ay3,dudsig_,dvdsig_;
	double *u,*v,*un,*vn,*ax,*ay;
    const int ID;
    
    // printing
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    double ddn;
    int fnpf_print_kinematicsprintcount;
    ofstream fout;

    // parallelisation
    double xstart,ystart,xend,yend;
};

#endif


