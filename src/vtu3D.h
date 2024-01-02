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
class force;
class vorticity;
class solver;
class probe_point;
class probe_pressure;
class probe_line;
class bedprobe_point;
class bedprobe_max;
class gage_discharge_x;
class gage_discharge_window_x;
class fsf_vtp;
class topo_vtp;
class cfd_state;
class bedshear_probe;
class bedshear_max;
class sloshing_force;
class print_porous;
class bedprobe_line_x;
class bedprobe_line_y;
class exportfile;
class flowfile_out;
class print_averaging;

#ifndef VTU3D_H_
#define VTU3D_H_

using namespace std;

class vtu3D : public printer, public nodefill 
{

public:
	vtu3D(lexer*,fdm*,ghostcell*);
	virtual ~vtu3D();
	virtual void start(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,multiphase*,sediment*);
    virtual void print_vtu(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,multiphase*,sediment*);
    virtual void print_stop(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,multiphase*,sediment*);
	virtual void ini(lexer*,fdm*,ghostcell*);

private:
    void print3D(fdm*,lexer*,ghostcell*,turbulence*,heat*,solver*,data*,concentration*,multiphase*,sediment*);
    void pvtu(fdm*,lexer*,ghostcell*,turbulence*,heat*,data*,concentration*,multiphase*,sediment*);
    void header(fdm*,lexer*,ghostcell*);
    void name_iter(fdm*,lexer*,ghostcell*);
    void name_time(fdm*,lexer*,ghostcell*);
    void piecename(fdm*,lexer*,ghostcell*, int);

    char name[200],pname[200],epsvar[200];
    int n,iin,offset[300];
    float ffn;
    int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double *printfsftime_wT;
    int *printfsfiter_wI;
    double phase;
    double zcoor;
	
	field5 eta;

    print_wsf *pwsf;
	print_wsf_theory *pwsf_theory;
    print_wsfline_x *pwsfline_x;
	print_wsfline_y *pwsfline_y;
    force **pforce;
    vorticity *pvort;
	probe_point *pprobe;
    probe_pressure *ppressprobe;
	probe_line *pline;
	bedprobe_point *pbedpt;
	bedprobe_line_x *pbedlinex;
	bedprobe_line_y *pbedliney;
	bedprobe_max *pbedmax;
	bedshear_probe *pbedshear;
	bedshear_max *pbedshearmax;
	gage_discharge_x *pq;
    gage_discharge_window_x *pqw;
	fsf_vtp *pfsf;
    topo_vtp *ptopo;
	cfd_state *pstate;
    sloshing_force *pslosh;
	print_porous *ppor;
    exportfile *pexport;
    flowfile_out *pflowfile;
    print_averaging *pmean;
};

#endif

