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
class fdm;
class ghostcell;
class slice;
class ioflow;

using namespace std;

#ifndef PROBE_VEL_THEORY_H_
#define PROBE_VEL_THEORY_H_

class probe_vel_theory : public boundarycheck
{
public:
    probe_vel_theory(lexer*,fdm*);
	virtual ~probe_vel_theory();

	void start(lexer*, fdm*, ghostcell*,ioflow*);


private:
    void ini_location(lexer*, fdm*);


	char name[100];

    int *iloc,*jloc,*kloc,*flag;
    int n,q;
	const int probenum;
    ofstream *pout;
	
	double uval,vval,wval,pval,kval,eval,edval;
    
    
    

};

#endif
