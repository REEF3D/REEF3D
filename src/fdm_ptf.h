/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/


#include"fdm.h"
#include"matrix2D.h"
#include"vec2D.h"
#include"cpt2D.h"

class lexer;

#ifndef FDM_PTF_H_
#define FDM_PTF_H_

using namespace std;

class fdm_ptf : public fdm
{
public:

    fdm_ptf(lexer*);
    ~fdm_ptf();

    // PTF
   // slice4 eta,eta_n,depth;
    slice4 Fifsf;
    slice4 K;
    slice4 WL;
    sliceint4 etaloc, breaking, breaklog, wet_n;
    slice4 coastline;
    slice4 vb;
    slice4 test2D;
    double *Fi_,*Uin_,*Uout_,*U_,*V_,*W_;
    double wd_criterion;
    
    matrix2D N;
    vec2D xvec,rvec;
    
    /*
    slice1 P;
    slice2 Q;
    
    slice4 bed;
    
	vec rhsvec;

	matrix_diag M;
	cpt C4,C4a,C6;

    double maxF,maxG,maxH;
    double wd_criterion;
	
	double t1,t2,t3,t4,t5;
     */
};

#endif
