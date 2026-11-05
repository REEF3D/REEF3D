/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef FDM_H_
#define FDM_H_

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
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"
#include"sliceint5.h"
#include"increment.h"
#include"vec.h"
#include"matrix_diag.h"
#include"looping.h"
#include<iostream>
#include<vector>

class lexer;

using namespace std;

// fields a grid of the mesh refinement leaves unallocated (cfd_amr::patch_objects): the patch
// objects never touch them (fdm(lexer*, unsigned lean))
enum : unsigned
{
    FDM_LEAN_FI     = 1u<<0,    // Fi (initialisation of the potential flow)
    FDM_LEAN_TEST   = 1u<<1,    // test (output)
    FDM_LEAN_CONC   = 1u<<2,    // conc (sediment suspension, CPM, output)
    FDM_LEAN_VOF    = 1u<<3,    // vof (VOF / PLIC)
    FDM_LEAN_NODE   = 1u<<4,    // nodeval (VTU output)
    FDM_LEAN_POR    = 1u<<5,    // porpart, porA, porB (VRANS)
    FDM_LEAN_WALLD  = 1u<<6     // walld (RANS)
};

class fdm : public increment
{
public:

    fdm(lexer*, unsigned lean=0);
    ~fdm();

	double gi,gj,gk;

	field1 u,F,Fext;
	field2 v,G,Gext;
	field3 w,H,Hext;
	field4 press;
    field4 Fi;
	field4 eddyv;
	field4 L;
	field4 ro,visc;
	field4 phi;
	field4 conc;
    field4 test;
	field4a topo,solid;
	field4a fb;
	field4a porosity,porpart,porA,porB;
	field5 walld;
	 
	fieldint5 nodeval;
    sliceint5 nodeval2D;
   
    // 6DOF
    field1 fbh1;
    field2 fbh2;
    field3 fbh3;
    field4 fbh4;
    field4 fbh5;
    
    //PLIC
    field4 nX,nY,nZ,Alpha;
    field4 phasemarker;
    field4 vof;
    field4 vof_nt,vof_nb,vof_st,vof_sb; //2D
    field4 vof_nte,vof_ntw,vof_nbe,vof_nbw,vof_ste,vof_stw,vof_sbe,vof_sbw; //3D
    
    // PTF
    slice4 eta,eta_n,depth,WL;
    slice4 Fifsf;
    slice4 K;
    sliceint4 etaloc;
    
    slice1 P;
    slice2 Q;
    
    slice4 bed;
    
	vec rhsvec;

	matrix_diag M;
    std::vector<unsigned long long> ckids;     // checkpoints (ckpt_registry.h)

    double maxF,maxG,maxH;
    double wd_criterion;
	
	double t1,t2,t3,t4,t5;
};

#endif
