/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"onephase.h"
#include"increment.h"
#include"slice4.h"

class field;
class fnpf_convection;

using namespace std;

#ifndef FNPF_FG_BED_UPDATE_H_
#define FNPF_FG_BED_UPDATE_H_

class fnpf_fg_bed_update : public increment
{
public:
    fnpf_fg_bed_update(lexer*, fdm*, ghostcell*);
	virtual ~fnpf_fg_bed_update();
    
	virtual void bedbc(lexer*, fdm*, ghostcell*,field&);
    virtual void waterdepth(lexer*, fdm*, ghostcell*);

private: 
    
    fnpf_convection *pconvec;


};

#endif
