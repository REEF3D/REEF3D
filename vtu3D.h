/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
class print_wsfline;
class print_wsfline_y;
class force;
class forcesolid;
class vorticity;
class solver;
class probe_point;
class probe_line;
class print_runup;
class bedprobe_point;
class bedprobe_max;
class gage_discharge;
class fsf_vtp;
class state;
class bedshear_probe;
class bedshear_max;
class sloshing_force;
class print_porous;
class bedprobe_line_x;
class bedprobe_line_y;
class exportfile;
class flowfile_out;

#ifndef VTU3D_H_
#define VTU3D_H_

using namespace std;

class vtu3D : public printer, public nodefill 
{

public:
	vtu3D(lexer*,fdm*,ghostcell*);
	virtual ~vtu3D();
	virtual void start(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,sediment*);
    virtual void print_vtu(fdm*,lexer*,ghostcell*,turbulence*,heat*,ioflow*,solver*,data*,concentration*,sediment*);
	virtual void ini(lexer*,fdm*,ghostcell*);

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
    print_wsfline *pwsfline;
	print_wsfline_y *pwsfline_y;

    force **pforce;
    forcesolid **pforcesolid;
    vorticity *pvort;
	probe_point *pprobe;
	probe_line *pline;
	print_runup *prunup;
	bedprobe_point *pbedpt;
	bedprobe_line_x *pbedlinex;
	bedprobe_line_y *pbedliney;
	bedprobe_max *pbedmax;
	bedshear_probe *pbedshear;
	bedshear_max *pbedshearmax;
	gage_discharge *pq;
	fsf_vtp *pfsf;
	state *pstate;
    sloshing_force *pslosh;
	print_porous *ppor;
    exportfile *pexport;
    flowfile_out *pflowfile;
};

#endif

