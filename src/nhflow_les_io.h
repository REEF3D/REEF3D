/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_turbulence.h"
#include"nhflow_strain.h"
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

#ifndef NHFLOW_LES_IO_H_
#define NHFLOW_LES_IO_H_

class nhflow_les_io : public nhflow_turbulence, public nhflow_strain
{
public:
    nhflow_les_io(lexer*,fdm_nhf*);
	virtual ~nhflow_les_io();

    virtual void print_3D(lexer*, fdm_nhf*, ghostcell*,ofstream&);
    virtual void ini(lexer*, fdm_nhf*, ghostcell*);
    virtual void plain_wallfunc(lexer*, fdm_nhf*, ghostcell*);
    virtual void inflow(lexer*, fdm_nhf*, ghostcell*);
    virtual double kinval(int,int,int);
    virtual double epsval(int,int,int);
	virtual void gcupdate(lexer*, fdm_nhf*, ghostcell*);
	virtual double ccipol_kinval(lexer*,ghostcell*,double,double,double);
	virtual double ccipol_epsval(lexer*,ghostcell*,double,double,double);
    virtual double ccipol_a_kinval(lexer*,ghostcell*,double,double,double);
	virtual double ccipol_a_epsval(lexer*,ghostcell*,double,double,double);
    virtual void kinget(int,int,int,double);
    virtual void epsget(int,int,int,double);
    
    virtual void isource(lexer*,fdm_nhf*);
	virtual void jsource(lexer*,fdm_nhf*);
	virtual void ksource(lexer*,fdm_nhf*);

    virtual void name_pvtu(lexer*, fdm_nhf*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &);
	

private:
    void tau_calc(fdm_nhf*, lexer*, double);
    void kepsini_default(lexer*,fdm_nhf*,ghostcell*);

	float ffn;
	int n,iin;
    int ii,jj,kk;
};

#endif


