/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
class fdm;
class ghostcell;
class onephase;
class ioflow;
class slice;
class field;

using namespace std;

#ifndef FNPF_FG_FSF_UPDATE_H_
#define FNPF_FG_FSF_UPDATE_H_

class fnpf_fg_fsf_update : public increment
{
public:
    fnpf_fg_fsf_update(lexer*, fdm*, ghostcell*);
	virtual ~fnpf_fg_fsf_update();
    
    virtual void fsfepol(lexer*, fdm*, ghostcell*,slice&,field&);
	virtual void fsfupdate(lexer*, fdm*, ghostcell*,ioflow*,onephase*,slice&);
    virtual void etaloc_sig(lexer*, fdm*, ghostcell*);
    virtual void etaloc(lexer*, fdm*, ghostcell*);
    virtual void fsfbc_sig(lexer*, fdm*, ghostcell*,slice&,field&);
    virtual void fsfbc(lexer*, fdm*, ghostcell*,slice&,field&);
    
    void velcalc(lexer*, fdm*, ghostcell *pgc, field&);
    void velcalc_sig(lexer*, fdm*, ghostcell *pgc, field&);
    
private: 
    int gcval,gcval_u,gcval_v,gcval_w;

};

#endif
