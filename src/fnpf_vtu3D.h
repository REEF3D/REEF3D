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

#include"printer.h"
#include"increment.h"

class fdm_fnpf;
class force_ale;
class fnpf_print_wsf;
class fnpf_print_wsf_theory;
class fnpf_print_wsfline;
class fnpf_print_wsfline_y;
class fnpf_vtp_fsf;
class fnpf_vtp_bed;
class fnpf_state;
class fnpf_breaking_log;
class fnpf_print_Hs;
class potentialfile_out;
class ioflow;

#ifndef FNPF_VTU3D_H_
#define FNPF_VTU3D_H_

using namespace std;

class fnpf_vtu3D : public increment
{

public:
	fnpf_vtu3D(lexer*,fdm_fnpf*,ghostcell*);
	virtual ~fnpf_vtu3D();
	virtual void start(lexer*,fdm_fnpf*,ghostcell*,ioflow*);
    virtual void print_stop(lexer*,fdm_fnpf*,ghostcell*);
    virtual void print_vtu(lexer*,fdm_fnpf*,ghostcell*);
    
private:
    void pvtu(lexer*,ghostcell*);
    void name_iter(lexer*,ghostcell*);
    void name_time(lexer*,ghostcell*);
    void piecename(lexer*,ghostcell*, int);

    char name[200],pname[200],epsvar[200];
    int n,iin,offset[200];
    float ffn;
    int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double *printfsftime_wT;
    int *printfsfiter_wI;
    double phase;
    double zcoor;
    
    int printcount;
    
    fnpf_print_wsf *pwsf;
    fnpf_print_wsf_theory *pwsf_theory;
    fnpf_print_wsfline *pwsfline;
    fnpf_print_wsfline_y *pwsfline_y;
    potentialfile_out *ppotentialfile;
    fnpf_vtp_fsf *pfsf;
    fnpf_vtp_bed *pbed;
    fnpf_state *pstate;
    fnpf_breaking_log *pbreaklog;
	force_ale **pforce_ale;
    fnpf_print_Hs *phs;
};

#endif

