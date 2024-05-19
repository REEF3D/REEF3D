/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class ioflow;

using namespace std;

#ifndef NHFLOW_VEL_PROBE_THEORY_H_
#define NHFLOW_VEL_PROBE_THEORY_H_

class nhflow_vel_probe_theory : public boundarycheck
{
public:
    nhflow_vel_probe_theory(lexer*,fdm_nhf*);
	virtual ~nhflow_vel_probe_theory();

	void start(lexer*, fdm_nhf*, ghostcell*,ioflow*);


private:
    void ini_location(lexer*, fdm_nhf*);


	char name[100];

    int *iloc,*jloc,*kloc,*flag;
    int n,q;
	const int probenum;
    ofstream *pout;
	
	double uval,vval,wval,pval,kval,eval,edval;
    
    
    

};

#endif
