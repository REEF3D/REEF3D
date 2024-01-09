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

#include"nhflow_turbulence.h"
#include"nhflow_strain.h"
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;

using namespace std;

#ifndef NHFLOW_RANS_IO_H_
#define NHFLOW_RANS_IO_H_

class nhflow_rans_io : public nhflow_turbulence, public nhflow_strain
{
public:
    nhflow_rans_io(lexer*,fdm_nhf*);
	virtual ~nhflow_rans_io();

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

    virtual void name_pvtu(lexer*, fdm_nhf*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &);
    
    double *KIN,*EPS,*EDDYV0;
    int *wallf;
	
	double const ke_c_1e, ke_c_2e,ke_sigma_k,ke_sigma_e;
	double const kw_alpha, kw_beta,kw_sigma_k,kw_sigma_w;
	double const sst_alpha1, sst_alpha2, sst_beta1, sst_beta2, sst_sigma_k1, sst_sigma_k2, sst_sigma_w1, sst_sigma_w2;

private:
    void tau_calc(fdm_nhf*, lexer*, double);
    void kepsini_default(lexer*,fdm_nhf*,ghostcell*);

	float ffn;
	int q,iin;
	int gcval_kin,gcval_eps,gcval_edv;

	double M,I,tau,H,B,ks,kinbed,uvel,refwalldist,fc;
	double kinw,epsw;
	double walld,ddn,depth;
};

#endif


