/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"fnpf_fg.h"
#include"slice4.h"

class fnpf_fg_laplace;
class field;
class fnpf_convection;

using namespace std;

#ifndef FNPF_FG_FSFBC_H_
#define FNPF_FG_FSFBC_H_

class fnpf_fg_fsfbc : public increment
{
public:
	fnpf_fg_fsfbc(lexer*, fdm*, ghostcell*);
	virtual ~fnpf_fg_fsfbc();
    
    
    void fsfdisc(lexer*,fdm*,ghostcell*,slice&,slice&,field&);
    void kfsfbc(lexer*,fdm*,ghostcell*);
    void dfsfbc(lexer*,fdm*,ghostcell*,slice&);
    

    fnpf_convection *pconvec;

    double ivel,jvel,kvel;
    
    slice4 Fx,Fy,Fz;
    slice4 Ex,Ey;

};

#endif
