/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef PTF_FSF_UPDATE_H_
#define PTF_FSF_UPDATE_H_

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class ioflow;
class slice;
class field;

using namespace std;

class ptf_fsf_update : public increment
{
public:
    ptf_fsf_update(lexer*, fdm*, ghostcell*);
	virtual ~ptf_fsf_update();
    
    void fsfepol(lexer*, fdm*, ghostcell*,slice&,field&);
	void fsfupdate(lexer*, fdm*, ghostcell*,ioflow*,slice&);
    void etaloc(lexer*, fdm*, ghostcell*);
    void fsfbc(lexer*, fdm*, ghostcell*,slice&,field&);
    void fsfbc0(lexer*, fdm*, ghostcell*,slice&,field&);
    void fsfbc1(lexer*, fdm*, ghostcell*,slice&,field&);
    
    void velcalc(lexer*, fdm*, ghostcell *pgc, field&);
    
private: 
    int gcval,gcval_u,gcval_v,gcval_w;

};

#endif
