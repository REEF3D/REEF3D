/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"fnpf_sg.h"
#include"slice4.h"

class fnpf_sg_laplace;
class field;
class fnpf_discrete;
class fnpf_ddx;

using namespace std;

#ifndef FNPF_SG_FSFBC_H_
#define FNPF_SG_FSFBC_H_

class fnpf_sg_fsfbc : public increment 
{
public:
	fnpf_sg_fsfbc(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_sg_fsfbc();
    
    void fsfdisc(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&);
    void kfsfbc(lexer*,fdm_fnpf*,ghostcell*);
    void dfsfbc(lexer*,fdm_fnpf*,ghostcell*,slice&);
    void fsfwvel(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&);

    fnpf_discrete *pdisc;
    fnpf_ddx *pddx;

    double ivel1,ivel2,jvel1,jvel2,kvel;
    
    slice4 Fx,Fy;
    slice4 Ex,Ey;
    slice4 Exx,Eyy;
    slice4 Bx,By;
    slice4 Bxx,Byy;
};

#endif
