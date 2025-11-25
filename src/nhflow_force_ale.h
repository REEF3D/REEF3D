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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef NHFLOW_FORCE_ALE_H_
#define NHFLOW_FORCE_ALE_H_

#include"nhflow_gradient.h"
#include"fieldint5.h"
#include"field5.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

class nhflow_force_ale :  public nhflow_gradient
{

public:
	nhflow_force_ale(lexer*,fdm_nhf*,ghostcell*,int);
	virtual ~nhflow_force_ale();
	void start(lexer*,fdm_nhf*,ghostcell*) override;
    void ini(lexer*,fdm_nhf*,ghostcell*) override;

private:	
	
    void force_ale_force(lexer*,fdm_nhf*,ghostcell*);
	void print_force_ale(lexer*,fdm_nhf*,ghostcell*);
    void print_ini(lexer*,fdm_nhf*,ghostcell*);
	double dndt_f(lexer*, fdm_nhf*, ghostcell*);
	double dudsig_f(lexer*, fdm_nhf*, ghostcell*);
	double dvdsig_f(lexer*, fdm_nhf*, ghostcell*);
	double dudxi(lexer*, fdm_nhf*, ghostcell*);
	double dvdxi(lexer*, fdm_nhf*, ghostcell*);
	
    // force ale variabes
    double Fx1,Fy1,Fx,Fy,xc,yc,rc,cd,cm,dtn,eta2n,ax1,ay1,ax2,ay2,ax3,ay3,ax,ay;
	double *un, *unn, *vn,*vnn;
    const int ID;
    
    double grad,dndt,dudsig,dvdsig;
    double simtime_n,simtime_nn;
    double eta_n,eta_nn;
    
    // printing
    char name[100],pname[100],epsvar[100];
    int iin,offset[100];
    float ffn;
    int force_aleprintcount;
    ofstream fout;

    // parallelisation
    double xstart,ystart,xend,yend;
};

#endif


