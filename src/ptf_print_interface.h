/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"ptf_printer.h"
#include"ptf_nodefill.h"
#include"field5.h"

class topo;
class ptf_print_wsf;
class ptf_print_wsf_theory;
class ptf_print_wsfline_x;
class ptf_print_wsfline_y;
class solver_ptf;
class ptf_probe_point;
class ptf_probe_line;
class ptf_fsf_vtp;
class ptf_state;

#ifndef PTF_PRINT_INTERFACE_H_
#define PTF_PRINT_INTERFACE_H_

using namespace std;

class print_interface : public printer, public nodefill 
{

public:
	print_interface(lexer*,fdm_ptf*,ghostcell*);
	virtual ~print_interface();
	virtual void start(fdm_ptf*,lexer*,ghostcell*,ioflow*,solver*);
	

private:
    void print3D(fdm_ptf*,lexer*,ghostcell*,solver*);
    void pvtu(fdm_ptf*,lexer*,ghostcell*);
    void header(fdm_ptf*,lexer*,ghostcell*);
    void name_iter(fdm_ptf*,lexer*,ghostcell*);
    void name_time(fdm_ptf*,lexer*,ghostcell*);
    void piecename(fdm_ptf*,lexer*,ghostcell*, int);
	void ggcfacet_fill(lexer*,fdm_ptf*,ghostcell*,field&);

    char name[200],pname[200],epsvar[200];
    int n,iin,offset[200];
    float ffn;
    int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double phase;
	
	field5 eta;

    ptf_print_wsf *pwsf;
	ptf_print_wsf_theory *pwsf_theory;
    ptf_print_wsfline_x *pwsfline_x;
	ptf_print_wsfline_y *pwsfline_y;

	ptf_probe_point *pprobe;
	ptf_probe_line *pline;
	ptf_fsf_vtp *pfsf;
	ptf_state *pstate;
    

#endif