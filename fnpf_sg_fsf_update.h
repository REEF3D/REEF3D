/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"increment.h"

class lexer;
class fdm_fnpf;
class ghostcell;
class onephase;
class ioflow;
class slice;
class field;
class vec;

using namespace std;

#ifndef FNPF_SG_FSF_UPDATE_H_
#define FNPF_SG_FSF_UPDATE_H_

class fnpf_sg_fsf_update : public increment
{
public:
    fnpf_sg_fsf_update(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_sg_fsf_update();
    
    virtual void fsfepol(lexer*, fdm_fnpf*, ghostcell*,slice&,field&);
	virtual void fsfupdate(lexer*, fdm_fnpf*, ghostcell*,ioflow*,onephase*,slice&);
    virtual void etaloc_sig(lexer*, fdm_fnpf*, ghostcell*);
    virtual void etaloc(lexer*, fdm_fnpf*, ghostcell*);
    virtual void fsfbc_sig(lexer*, fdm_fnpf*, ghostcell*,slice&,double*);
    virtual void fsfbc(lexer*, fdm_fnpf*, ghostcell*,slice&,field&);
    
    void velcalc(lexer*, fdm_fnpf*, ghostcell *pgc, field&);
    void velcalc_sig(lexer*, fdm_fnpf*, ghostcell *pgc, double*);
    
private: 
    int gcval,gcval_u,gcval_v,gcval_w;


};

#endif
