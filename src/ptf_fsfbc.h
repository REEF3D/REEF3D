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
class ptf_coastline;

using namespace std;

#ifndef PTF_FSFBC_H_
#define PTF_FSFBC_H_

class ptf_fsfbc : public increment
{
public:
	ptf_fsfbc(lexer*, fdm_ptf*, ghostcell*);
	virtual ~ptf_fsfbc();
    
    
    void fsfdisc(lexer*,fdm_ptf*,ghostcell*,slice&,slice&,field&);
    void fsfdisc_ini(lexer*,fdm_ptf*,ghostcell*,slice&,slice&);
    void kfsfbc(lexer*,fdm_ptf*,ghostcell*);
    void dfsfbc(lexer*,fdm_ptf*,ghostcell*,slice&);
    void fsfwvel(lexer*,fdm_ptf*,ghostcell*,slice&,slice&);
    double fz(lexer*,fdm_ptf*,field&,slice&);
    
    void breaking(lexer*,fdm_ptf*,ghostcell*,slice&,slice&,slice&,double);
    void breaking_wd(lexer*,fdm_ptf*,ghostcell*,slice&,slice&,slice&,double);
    void damping(lexer*,fdm_ptf*,ghostcell*,slice&,int,double);
    void damping_wd(lexer*,fdm_ptf*,ghostcell*,slice&,int,double);
    void filter(lexer*, fdm_ptf*,ghostcell*, slice&);
    void filter_wd(lexer*, fdm_ptf*,ghostcell*, slice&);
    void coastline_eta(lexer*,fdm_ptf*,ghostcell*,slice&);
    void coastline_fi(lexer*,fdm_ptf*,ghostcell*,slice&);
    void wetdry(lexer*,fdm_ptf*,ghostcell*,slice&,slice&);

    fnpf_convection *pconvec;
    fnpf_convection *pdx;
    ptf_coastline *pcoast;

    double ivel,jvel,kvel;
    
    slice4 Fx,Fy,Fz;
    slice4 Ex,Ey;
    slice4 Bx,By;
    solver2D *psolv;
    double visc;
    double grad, teta;
    
private:
    double rb3(lexer*,double);
    double rb4(lexer*,double);
    double dist3,dist4,expinverse,db;
    sliceint4 bx,by;
    int count_n;

};

#endif
