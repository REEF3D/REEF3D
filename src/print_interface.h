/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"printer.h"
#include"nodefill.h"
#include"field5.h"

class turbulence;
class heat;
class suspended;
class bedload;
class topo;
class print_wsf;
class print_wsf_theory;
class print_wsfline_x;
class print_wsfline_y;
class force;
class vorticity;
class solver;
class probe_point;
class probe_line;
class bedprobe_point;
class bedprobe_max;
class gage_discharge_x;
class fsf_vtp;
class state;
class bedshear_probe;
class sloshing_force;
class print_porous;

#ifndef PRINT_INTERFACE_H_
#define PRINT_INTERFACE_H_

using namespace std;

class print_interface : public printer, public nodefill 
{

public:
	print_interface(lexer*,fdm*,ghostcell*);
	virtual ~print_interface();
	virtual void start(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,sediment*);
	

private:
    void print3D(fdm*,lexer*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,sediment*);
    void pvtu(fdm*,lexer*,ghostcell*,turbulence*,heat*,data*,concentration*,sediment*);
    void header(fdm*,lexer*,ghostcell*);
    void name_iter(fdm*,lexer*,ghostcell*);
    void name_time(fdm*,lexer*,ghostcell*);
    void piecename(fdm*,lexer*,ghostcell*, int);
	void ggcfacet_fill(lexer*,fdm*,ghostcell*,field&);

    char name[200],pname[200],epsvar[200];
    int n,iin,offset[200];
    float ffn;
    int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double phase;
	
	field5 eta;

    print_wsf *pwsf;
	print_wsf_theory *pwsf_theory;
    print_wsfline_x *pwsfline_x;
	print_wsfline_y *pwsfline_y;

    vorticity *pvort;
	probe_point *pprobe;
	probe_line *pline;
	bedprobe_point *pbedpt;
	bedprobe_max *pbedmax;
	bedshear_probe *pbedshear;
	gage_discharge_x *pq;
	fsf_vtp *pfsf;
	state *pstate;
    sloshing_force *pslosh;
	print_porous *ppor;
};

#endif

