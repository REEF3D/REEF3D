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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"pressure.h"
#include"increment.h"
#include"field4.h"

class heat;
class concentration;
class density;

using namespace std;

#ifndef PJM_IMEX_H_
#define PJM_IMEX_H_

class pjm_IMEX : public pressure, public increment
{

public:

	pjm_IMEX(lexer*, fdm*, heat*&, concentration*&);
	virtual ~pjm_IMEX();

	virtual void start(fdm*,lexer* , poisson*, solver*, ghostcell*, ioflow*, field&, field&, field&,double);
	virtual void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void upgrad(lexer*,fdm*,slice&,slice&);
	virtual void vpgrad(lexer*,fdm*,slice&,slice&);
	virtual void wpgrad(lexer*,fdm*,slice&,slice&);
	virtual void ucorr(lexer*,fdm*,field&,double);
	virtual void vcorr(lexer*,fdm*,field&,double);
	virtual void wcorr(lexer*,fdm*,field&,double);
    
    field4 pcorr, Fp, pressn;

private:

	void debug(lexer*,fdm*);
	void pressure_norm(lexer*,fdm*,ghostcell*);
    void presscorr(lexer*p,fdm *a,ghostcell*,field&,field&,field&,field&, double);
	void vel_setup(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	
    double starttime, endtime, sum_old, omega_k;
	int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
    
    density *pd;
};



#endif

