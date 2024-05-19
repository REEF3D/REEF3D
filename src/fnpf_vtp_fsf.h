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
class fdm_fnpf;
class ghostcell;
class ioflow;

using namespace std;

#ifndef FNPF_VTP_FSF_H_
#define FNPF_VTP_FSF_H_

class fnpf_vtp_fsf : public increment
{
public:
	fnpf_vtp_fsf(lexer*,fdm_fnpf*,ghostcell*);
	virtual ~fnpf_vtp_fsf();
	
    virtual void start(lexer*,fdm_fnpf*,ghostcell*);
    virtual void print2D(lexer*,fdm_fnpf*,ghostcell*);
	
private:
	
	void etend(lexer*,fdm_fnpf*,ghostcell*);
	void pvtu(lexer*,fdm_fnpf*,ghostcell*);
	void name_iter(lexer*,fdm_fnpf*,ghostcell*);
    void piecename(lexer*,fdm_fnpf*,ghostcell*,int);
	
	
	char name[200],pname[200];
    int n,iin,offset[200];
    float ffn;
	
	double xs_local,ys_local,zs_local,xe_local,ye_local,ze_local;
	double xs_global,ys_global,zs_global,xe_global,ye_global,ze_global;
    
    int gcval_eta, gcval_fifsf;
    int printcount;
	

};

#endif
