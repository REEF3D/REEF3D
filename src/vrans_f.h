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

#ifndef VRANS_F_H_
#define VRANS_F_H_

#include"vrans.h"
#include"increment.h"
#include"field4a.h"

using namespace std;

class vrans_f : public vrans, public increment
{
public:
	vrans_f(lexer*, ghostcell*);
	virtual ~vrans_f();

	void initialize_cfd(lexer*, fdm*, ghostcell*) override final;	
	void start(lexer*, fdm*, ghostcell*, int) override final {};
    void sed_update(lexer*, fdm*, ghostcell*) override final;	
    void sedpart_update(lexer*, fdm*, ghostcell*, field&, field&) override final;
	
	void u_source(lexer*, fdm*) override final;
	void v_source(lexer*, fdm*) override final;
	void w_source(lexer*, fdm*) override final;
    
    void ke_source(lexer*, fdm*, field&) override final;
    void kw_source(lexer*, fdm*, field&) override final;
    void eps_source(lexer*, fdm*, field&, field&) override final;
    void omega_source(lexer*, fdm*, field&, field&) override final;
    
    void eddyv_func(lexer*, fdm*) override final;
    
    void veltimesave(lexer*,fdm*,ghostcell*) override final;
	
private:
	
	field4a alpha,beta;
	
	double Apor(double,double,double,double);
	double Bpor(double,double,double);
	
	int count;
    
    double Aporval,Bporval,porval,partval,alphaval,betaval,viscval;
	double val;
	double porousterm;
	const double Cval;
};

#endif
