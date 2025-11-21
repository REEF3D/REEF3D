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

#ifndef SFLOW_VTP_FSF_H_
#define SFLOW_VTP_FSF_H_

#include"increment.h"
#include"vtp3D.h"

class lexer;
class fdm2D;
class ghostcell;
class ioflow;
class sediment;
class sflow_print_wsf;
class sflow_print_wsf_theory;
class sflow_print_wsfline;
class sflow_print_wsfline_y;
class sflow_print_probe_da;
class sflow_turbulence;
class sflow_state;
class fnpf_print_Hs;

using namespace std;

class sflow_vtp_fsf : public increment, private vtp3D
{
public:
    sflow_vtp_fsf(lexer*,fdm2D*,ghostcell*);
    virtual ~sflow_vtp_fsf() = default;

    void start(lexer*,fdm2D*,ghostcell*,ioflow*,sflow_turbulence*,sediment*);
    void print2D(lexer*,fdm2D*,ghostcell*,sflow_turbulence*,sediment*);

private:
    void pvtp(lexer*,fdm2D*,ghostcell*,sflow_turbulence*,sediment*,int);

    char name[200];
    int n,iin,offset[200];
    float ffn;

    sflow_print_wsf *pwsf;
    sflow_print_wsf_theory *pwsf_theory;
    sflow_print_wsfline *pwsfline;
    sflow_print_wsfline_y *pwsfline_y;
    sflow_print_probe_da *pprobe;
    sflow_state *pstate;
    sflow_state *pstate_restart;
    fnpf_print_Hs *phs;
};

#endif
