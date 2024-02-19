/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"ptf.h"
#include"slice4.h"
#include"sliceint4.h"

class ptf_laplace;
class field;
class fnpf_convection;
class solver2D;

using namespace std;

#ifndef PTF_FSFBC_H_
#define PTF_FSFBC_H_

class ptf_fsfbc : public increment
{
public:
	ptf_fsfbc(lexer*, fdm_ptf*, ghostcell*);
	virtual ~ptf_fsfbc();
    
    
    void fsfdisc(lexer*,fdm_ptf*,ghostcell*,slice&,slice&,field&);
    void kfsfbc(lexer*,fdm_ptf*,ghostcell*);
    void dfsfbc(lexer*,fdm_ptf*,ghostcell*,slice&);
    void fsfwvel(lexer*,fdm_ptf*,ghostcell*,slice&,slice&);
    double fz(lexer*,fdm_ptf*,field&,slice&);
    
    void breaking(lexer*,fdm_ptf*,ghostcell*,slice&,slice&,slice&,double);
    virtual void damping(lexer*,fdm_ptf*,ghostcell*,slice&,int,double);
    void filter(lexer*, fdm_ptf*,ghostcell*, slice&);

    fnpf_convection *pconvec;

    double ivel,jvel,kvel;
    
    slice4 Fx,Fy,Fz;
    slice4 Ex,Ey;
    sliceint4 bx,by;
    solver2D *psolv;
    double visc;
    double grad, teta;

};

#endif
