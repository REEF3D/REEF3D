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

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm;
class ghostcell;
class field;
class turbulence;

using namespace std;

#ifndef PROBE_LINE_H_
#define PROBE_LINE_H_

class probe_line : public boundarycheck
{
public:
    probe_line(lexer*,fdm*,ghostcell*);
	virtual ~probe_line();

	void start(lexer*, fdm*, ghostcell*,turbulence*);


private:
	void ini_global_location(lexer*, fdm*, ghostcell*);
    void ini_location(lexer*, fdm*, ghostcell*);
    void write(lexer*, fdm*, ghostcell*);
    int conv(double);
	char name[100];

    int **flag,**active,**active_all,**flag_all,**elid,**elid_all;
    int **iloc, **iloc_all, **kloc, **kloc_all;
    double **xs,**xs_all;
    
	int *totelnum,*totelnum_all,*elnum,**elnum_all,**displ;
	int maxelnum,maxlocalelnum,count,allcount;
	double *length,*ds,*norm;
	double **U,**V,**W,**P,**K,**E,**VT,**LS;
	double **U_all,**V_all,**W_all,**P_all,**K_all,**E_all,**VT_all,**LS_all,**VAL;

    int n,q;
	const int probenum;
	const double eps;
    ofstream *lineout;
	
	double domain_xs,domain_xe,domain_ys,domain_ye,domain_zs,domain_ze;
	
	double uval,vval,wval,pval,kval,eval,edval;
	int check;
	double t;
	double xloc,yloc,zloc;
	double xp,yp,zp;
	double eps_xs, eps_xe, eps_ys, eps_ye, eps_zs, eps_ze;
	int linecount;

};

#endif
