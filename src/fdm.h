/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"field4a.h"
#include"field5.h"
#include"fieldint5.h"
#include"fieldint1.h"
#include"fieldint2.h"
#include"fieldint3.h"
#include"fieldint4.h"
#include"fieldint4a.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"
#include"sliceint5.h"
#include"increment.h"
#include"vec.h"
#include"matrix_diag.h"
#include"cpt.h"
#include"looping.h"
#include"iterators.h"
#include<iostream>
#include<vector>

class lexer;

#ifndef FDM_H_
#define FDM_H_

using namespace std;

class fdm : public increment
{
public:

    fdm(lexer*);

	double gi,gj,gk;

	field1 u,F;
	field2 v,G;
	field3 w,H;
	field4 press;
    field4 Fi;
	field4 eddyv;
	field4 L;
	field4 ro,visc;
	field4 phi,vof;
	field4 conc;
    field4 test;
	field4a topo,solid;
	field4a fb;
	field4a porosity;
	field5 walld;
	 
	fieldint5 nodeval;
    sliceint5 nodeval2D;
   
    // 6DOF
    field1 fbh1;
    field2 fbh2;
    field3 fbh3;
    field4 fbh4;
    field4 fbh5;
    
    
    // PTF
    slice4 eta,eta_n,depth;
    slice4 Fifsf;
    slice4 K;
    sliceint4 etaloc;
    
    slice1 P;
    slice2 Q;
    
    slice4 bed;
    
	vec rhsvec;

	matrix_diag M;
	cpt C4,C4a,C6;

    double maxF,maxG,maxH;
    double wd_criterion;
	
	double t1,t2,t3,t4,t5;
};

#endif
