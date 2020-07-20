/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
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

#ifndef PTF_FSF_UPDATE_H_
#define PTF_FSF_UPDATE_H_

class ptf_fsf_update : public increment
{
public:
    ptf_fsf_update(lexer*, fdm*, ghostcell*);
	virtual ~ptf_fsf_update();
    
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
