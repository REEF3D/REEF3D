/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf.h"
#include"increment.h"
#include"fnpf_fsf_update.h"
#include"fnpf_bed_update.h"
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;
class print;
class ioflow;
class reini;
class onephase;

using namespace std;

#ifndef FNPF_INI_H_
#define FNPF_INI_H_

class fnpf_ini : public fnpf, public increment, public fnpf_fsf_update, public fnpf_bed_update
{
public:
	fnpf_ini(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_ini();
    
    virtual void ini(lexer*, fdm_fnpf*, ghostcell*, ioflow*, reini*, onephase*);
    
    void velcalc(lexer*, fdm_fnpf*, ghostcell *pgc, field&);
    
    void fnpf_restart(lexer*, fdm_fnpf*, ghostcell *pgc);
    void fnpf_restart_mainheader(lexer*, fdm_fnpf*, ghostcell *pgc);
    void fnpf_restart_read(lexer*, fdm_fnpf*, ghostcell *pgc);
    void filename(lexer*, fdm_fnpf*, ghostcell *pgc,int);
    
private:

    void lsm_ini(lexer*, fdm_fnpf*, ghostcell*, ioflow*);

    int gcval,gcval_u,gcval_v,gcval_w;
    
    char name[500];
    
    int iin,file_type;
    float ffn;
	double ddn;
	int printcount;
    ifstream result;

};

#endif
