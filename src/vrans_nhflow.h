/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef VRANS_F_H_
#define VRANS_F_H_

#include"vrans.h"
#include"increment.h"

using namespace std;

class vrans_nhflow : public vrans, public increment
{
public:
	vrans_nhflow(lexer*, ghostcell*);
	virtual ~vrans_nhflow();

	void initialize(lexer*, fdm_nhf*, ghostcell*) override;	
	void start(lexer*, fdm_nhf*, ghostcell*, int) override {};
    void sed_update(lexer*, fdm_nhf*, ghostcell*) override;	
    void sedpart_update(lexer*, fdm_nhf*, ghostcell*, field&, field&) override;
	
	void u_source(lexer*, fdm_nhf*) override;
	void v_source(lexer*, fdm_nhf*) override;
	void w_source(lexer*, fdm_nhf*) override;
    
    void ke_source(lexer*, fdm_nhf*, field&) override;
    void kw_source(lexer*, fdm_nhf*, field&) override;
    void eps_source(lexer*, fdm_nhf*, field&, field&) override;
    void omega_source(lexer*, fdm_nhf*, field&, field&) override;
    
    void eddyv_func(lexer*, fdm_nhf*) override;
    
	
private:
	
	double *NPOR,*DPOR,*APOR,*BPOR;
	
	double Apor(double,double,double,double);
	double Bpor(double,double,double);
	
	int count;
    
    double Aporval,Bporval,porval,partval,alphaval,betaval,viscval;
	double val;
	double porousterm;
	const double Cval;
};

#endif
