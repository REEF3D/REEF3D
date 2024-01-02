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

#include"increment.h"

class lexer;
class fdm2D;
class ghostcell;
class ioflow;
class sediment;
class sflow_print_wsf;
class sflow_print_wsf_theory;
class sflow_print_wsfline;
class sflow_print_wsfline_y;
class sflow_print_bed;
class sflow_print_bedline;
class sflow_print_bedline_y;
class sflow_print_probe_da;
class sflow_turbulence;
class sflow_state;
class fnpf_print_Hs;

using namespace std;

#ifndef SFLOW_VTP_FSF_H_
#define SFLOW_VTP_FSF_H_

class sflow_vtp_fsf : public increment
{
public:
	sflow_vtp_fsf(lexer*,fdm2D*,ghostcell*);
	virtual ~sflow_vtp_fsf();

    virtual void start(lexer*,fdm2D*,ghostcell*,ioflow*,sflow_turbulence*,sediment*);
    virtual void print2D(lexer*,fdm2D*,ghostcell*,sflow_turbulence*,sediment*);

private:

	void etend(lexer*,fdm2D*,ghostcell*);
	void pvtp(lexer*,fdm2D*,ghostcell*,sflow_turbulence*,sediment*);
	void name_iter(lexer*,fdm2D*,ghostcell*);
    void piecename(lexer*,fdm2D*,ghostcell*,int);


	char name[200],pname[200];
    int n,iin,offset[200];
    float ffn;
    double ddn;

	double xs_local,ys_local,zs_local,xe_local,ye_local,ze_local;
	double xs_global,ys_global,zs_global,xe_global,ye_global,ze_global;

	sflow_print_wsf *pwsf;
	sflow_print_wsf_theory *pwsf_theory;
    sflow_print_wsfline *pwsfline;
    sflow_print_wsfline_y *pwsfline_y;
    sflow_print_probe_da *pprobe;
    sflow_print_bed *pbed;
    sflow_print_bedline *pbedline;
    sflow_print_bedline_y *pbedline_y;
    sflow_state *pstate;
    fnpf_print_Hs *phs;

};

#endif
