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
Author: Dave Kelly
--------------------------------------------------------------------*/

#include"increment.h"
#include"slice4.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice;

using namespace std;

#ifndef FNPF_PRINT_HS_H_
#define FNPF_PRINT_HS_H_

class fnpf_print_Hs : public increment
{
public:
    fnpf_print_Hs(lexer*,fdm_fnpf*);
	virtual ~fnpf_print_Hs();

	void start(lexer*, fdm_fnpf*, ghostcell*,slice&);
    
    slice4 ETAsum, ETAmean; //DKAF
    slice4 ETA2sum, ETAvar; //DKAF

private:
    int NumDT1;      
    double T_INTV_mean; // Averaging time for sig wave height
    double T_sum,dT_sum; 
    int wfcall;      
    double wtime;
    double stime;        // Start avreging after transients
    

};

#endif
