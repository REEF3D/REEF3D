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

#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;

using namespace std;

#ifndef NHFLOW_PRINT_WSF_H_
#define NHFLOW_PRINT_WSF_H_

class nhflow_print_wsf : public increment
{
public:
    nhflow_print_wsf(lexer*,fdm_nhf*);
	virtual ~nhflow_print_wsf();

	void height_gauge(lexer*, fdm_nhf*, ghostcell*,slice&);


private:
    void ini_location(lexer*, fdm_nhf*);
    void fill_eta(lexer*, fdm_nhf*, ghostcell*,slice&);
    void fill_deta(lexer*, fdm_nhf*, ghostcell*,slice&);
    void fill_Uhorz(lexer*, fdm_nhf*, ghostcell*,slice&);

	
	double *x,*y;
	int gauge_num;

    int *iloc,*jloc,*flag;
    double *wsf,*deta,*Uhorz;
    int n;
    ofstream wsfout,detaout,Uhorzout;
    

};

#endif
