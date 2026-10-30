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

#ifndef NHFLOW_TURBULENCE_H_
#define NHFLOW_TURBULENCE_H_

class fdm_nhf;
class lexer;
class nhflow_scalar_convection;
class nhflow_diffusion;
class solver;
class ghostcell;
class ioflow;
class vrans;
class vrans_nhflow;

#include<fstream>
#include<sstream>
#include<vector>
#include<functional>

using namespace std;

class nhflow_turbulence
{

public:
    // checkpoints (P 48, I 45): the arrays of the model (KIN, EPS, ...): fd(name, kind, double*), fi(.., int*)
    virtual void ckpt_fields(const std::function<void(const char*, char, double*)>&, const std::function<void(const char*, char, int*)>&) {}
    virtual ~nhflow_turbulence() {}

	virtual void start(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, vrans_nhflow*)=0;
    // the step in parts (k-epsilon, k-omega: 0 k, 1 epsilon/omega, 2 eddy viscosity); other models make the whole
    // step in part 0
    virtual void start_part(lexer *p, fdm_nhf *d, ghostcell *pgc, nhflow_scalar_convection *c, nhflow_diffusion *df, solver *s, ioflow *f, vrans_nhflow *v, int part)
    {
        if(part==0)
        start(p,d,pgc,c,df,s,f,v);
    }
    // set by nhflow_amr on a multi-block grid: takes over the step of level 0 (true) and runs the parts of all blocks
    std::function<bool(lexer*, fdm_nhf*, ghostcell*, nhflow_scalar_convection*, nhflow_diffusion*, solver*, ioflow*, vrans_nhflow*)> defer;
	virtual void ktimesave(lexer*, fdm_nhf*, ghostcell*)=0;
	virtual void etimesave(lexer*, fdm_nhf*, ghostcell*)=0;
	virtual void isource(lexer*, fdm_nhf*)=0;
	virtual void jsource(lexer*, fdm_nhf*)=0;
	virtual void ksource(lexer*,fdm_nhf*)=0;

	virtual void print_2D(lexer*, fdm_nhf*, ghostcell*,ofstream&,int)=0;
    virtual void print_3D(lexer*, fdm_nhf*, ghostcell*,std::vector<char>&, size_t&)=0;
    virtual void ini(lexer*, fdm_nhf*, ghostcell*)=0;
    virtual double kinval(int,int,int)=0;
    virtual double epsval(int,int,int)=0;
	virtual void gcupdate(lexer*, fdm_nhf*, ghostcell*)=0;
	virtual double ccipol_kinval(lexer*,ghostcell*,double,double,double)=0;
	virtual double ccipol_epsval(lexer*,ghostcell*,double,double,double)=0;
    virtual double ccipol_a_kinval(lexer*,ghostcell*,double,double,double)=0;
	virtual double ccipol_a_epsval(lexer*,ghostcell*,double,double,double)=0;
    virtual void kinget(int,int,int,double)=0;
    virtual void epsget(int,int,int,double)=0;

    virtual void name_ParaView_parallel(lexer*, ofstream&)=0;
    virtual void name_ParaView(lexer*, stringstream&, int*, int &)=0;
    virtual void offset_ParaView(lexer*, int*, int &)=0;

    virtual void name_pvtp(lexer*, fdm_nhf*, ghostcell*,ofstream&)=0;
    virtual void name_vtp(lexer*, fdm_nhf*, ghostcell*,ofstream&, int*, int &)=0;
    virtual void offset_ParaView_2D(lexer*, int*, int &)=0;

    double uref;
};

#endif
