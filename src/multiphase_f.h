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

class fdm;
class lexer;
class convection;
class solver;
class ghostcell;
class ioflow;
class reini;
class printer;
class field;
class freesurface;
class multiphase_fluid_update;
class heat;
class print_wsf;
class concentration;

#include"multiphase.h"
#include"field4.h"
#include<fstream>

using namespace std;

#ifndef MULTIPHASE_F_H_
#define MULTIPHASE_F_H_

class multiphase_f : public multiphase, public increment
{
public:
	multiphase_f(lexer*, fdm*, ghostcell*);
	virtual ~multiphase_f();
	virtual void start(lexer*,fdm*,ghostcell*,convection*,solver*,ioflow*,reini*,particle_corr*,printer*);
	virtual void ini(lexer*,fdm*,ghostcell*,ioflow*,printer*,convection*,solver*);
	virtual void update(lexer*,fdm*,ghostcell*);
	
	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void print_file(lexer*, fdm*, ghostcell*);
	virtual void nodefill(lexer*,fdm*,ghostcell*,field&);
    virtual double ls1val(int,int,int);
    virtual double ls2val(int,int,int);
	virtual double ccipol_ls1val(lexer*,ghostcell*,double,double,double);
	virtual double ccipol_ls2val(lexer*,ghostcell*,double,double,double);
    virtual void ls1get(int,int,int,double);
    virtual void ls2get(int,int,int,double);

    virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
	
	void logic(lexer*,fdm*,ghostcell*);
	
	freesurface *pfsf1,*pfsf2;
	reini *preini;
	multiphase_fluid_update *pupdate;
	heat *pheat;
	concentration *pconc;
	particle_corr *ppls;
	print_wsf *pwsf1;
	print_wsf *pwsf2;
	
	field4 ls1,ls2;
	
	int n;
	int iin;
	float ffn;
	double ddn;

private:
	double fx(double,double,double,double,double);
	double fz(double,double,double,double,double);
};

#endif