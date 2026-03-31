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

#ifndef VRANS_NHFLOW_F_H_
#define VRANS_NHFLOW_F_H_

#include"vrans_nhflow.h"
#include"nhflow_geometry.h"

using namespace std;

class vrans_nhflow_f final : public vrans_nhflow, public nhflow_geometry
{
public:
	vrans_nhflow_f(lexer*, fdm_nhf*, ghostcell*);
	virtual ~vrans_nhflow_f();

	void initialize(lexer*, fdm_nhf*, ghostcell*) override final;	
	void update(lexer*, fdm_nhf*, ghostcell*, int) override final;
	
	void u_source(lexer*, fdm_nhf*, slice&) override final;
	void v_source(lexer*, fdm_nhf*, slice&) override final;
	void w_source(lexer*, fdm_nhf*, slice&) override final;
    
    void ke_source(lexer*, fdm_nhf*, field&) override final;
    void kw_source(lexer*, fdm_nhf*, field&) override final;
    void eps_source(lexer*, fdm_nhf*, field&, field&) override final;
    void omega_source(lexer*, fdm_nhf*, field&, field&) override final;
    
    void eddyv_func(lexer*, fdm_nhf*) override final;
    
	
private:
    
    double Hporface(lexer*, fdm_nhf*, int, int, int);
	
	double *APOR,*BPOR;
	
	double Apor(double,double,double,double);
	double Bpor(double,double,double);
	
	int count;
    
    double H;
    
    double Aporval,Bporval,porval,partval,alphaval,betaval,viscval;
	double val;
	double porousterm;
	const double Cval;
};

#endif
