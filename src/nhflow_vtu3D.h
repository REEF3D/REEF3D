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

#ifndef NHFLOW_VTU3D_H_
#define NHFLOW_VTU3D_H_

#include"nhflow_printer.h"
#include"increment.h"

class fdm_nhf;
class ioflow;
class sediment;
class nhflow_print_wsf;
class nhflow_print_wsf_theory;
class nhflow_print_wsfline;
class nhflow_print_wsfline_y;
class nhflow_print_runup_gage_x;
class nhflow_print_runup_max_gage_x;
class nhflow_u_profile;
class nhflow_vtp_fsf;
class nhflow_vtp_bed;
class nhflow_state;
class nhflow_breaking_log;
class nhflow_vel_probe;
class nhflow_vel_probe_theory;
class nhflow_print_Hs;
class nhflow_turbulence;
class nhflow_force;
class nhflow_force_ale;
class bedshear_probe;
class bedshear_max;

using namespace std;

class nhflow_vtu3D : public nhflow_printer, public increment
{

public:
	nhflow_vtu3D(lexer*,fdm_nhf*,ghostcell*);
	virtual ~nhflow_vtu3D();
	virtual void start(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*);
    virtual void print_vtu(lexer*,fdm_nhf*,ghostcell*,nhflow_turbulence*,sediment*);
    virtual void print_stop(lexer*,fdm_nhf*,ghostcell*,ioflow*,nhflow_turbulence*,sediment*);
    
private:
    void pvtu(lexer*,fdm_nhf*,ghostcell*,nhflow_turbulence*,sediment*);
    void name_iter(lexer*,ghostcell*);
    void name_time(lexer*,ghostcell*);
    void piecename(lexer*,ghostcell*, int);

    char name[200],pname[200],epsvar[200];
    int n,iin,offset[200];
    float ffn;
    int jj;
    int gcval_phi,gcval_phiext;
	double *printtime_wT;
    double *printfsftime_wT;
    int *printfsfiter_wI;
    double phase;
    double zcoor;
    
    int printcount;
    
    nhflow_print_wsf *pwsf;
    nhflow_print_wsf_theory *pwsf_theory;
    nhflow_print_wsfline *pwsfline;
    nhflow_print_wsfline_y *pwsfline_y;
    nhflow_print_runup_gage_x *prunupx;
    nhflow_print_runup_max_gage_x *prunupmaxx;
    nhflow_u_profile *puprofile;
    nhflow_vtp_fsf *pfsf;
    nhflow_vtp_bed *pbed;
    nhflow_state *pstate;
    nhflow_breaking_log *pbreaklog;
    nhflow_vel_probe *pvel;
    nhflow_vel_probe_theory *pveltheo;
    nhflow_print_Hs *phs;
    nhflow_force **pforce;
    nhflow_force_ale **pforce_ale;
    
};

#endif

