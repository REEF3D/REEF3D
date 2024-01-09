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

#include"vrans.h"
#include"increment.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4a.h"

using namespace std;

#ifndef VRANS_VEG_H_
#define VRANS_VEG_H_

class vrans_veg : public vrans, public increment
{
public:
	vrans_veg(lexer*, ghostcell*);
	virtual ~vrans_veg();

	virtual void initialize(lexer*, fdm*, ghostcell*);	
	virtual void start(lexer*, fdm*, ghostcell*, net*&, int){};
    virtual void sed_update(lexer*, fdm*, ghostcell*);
	
	virtual void u_source(lexer*, fdm*);
	virtual void v_source(lexer*, fdm*);
	virtual void w_source(lexer*, fdm*);
    
    virtual void ke_source(lexer*, fdm*, field&);
    virtual void kw_source(lexer*, fdm*, field&);
    virtual void eps_source(lexer*, fdm*, field&, field&);
    virtual void omega_source(lexer*, fdm*, field&, field&);
    
    virtual void eddyv_func(lexer*, fdm*);
    
    virtual void veltimesave(lexer*,fdm*,ghostcell*);
	
private:
	
	field4a N,D,Cd;
    field1 un;
    field2 vn;
    field3 wn;
	
	double Apor(double,double,double,double);
	double Bpor(double,double,double);
	
	int count;
    
    double N_val,D_val,Cd_val,Fi,Fd;
    double porval,partval,alphaval,betaval,viscval;
	double val;
	double porousterm;
	const double Cval;
};

#endif
