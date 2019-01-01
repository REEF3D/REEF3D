/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"pressure.h"
#include"density.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

using namespace std;

#ifndef SIMPLE_H_
#define SIMPLE_H_

class simple : public pressure,  public density
{

public:
	simple(lexer* p, fdm *a);
	virtual ~simple();

	virtual void start(fdm*,lexer* p, poisson*, solver*, ghostcell*,momentum*,ioflow*, field&, field&, field&,double);
	virtual void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void vel_setup(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void ucorr(fdm*,lexer*p,field&);
	virtual void vcorr(fdm*,lexer*p,field&);
	virtual void wcorr(fdm*,lexer*p,field&);
	virtual void upgrad(lexer*,fdm*);
	virtual void vpgrad(lexer*,fdm*);
	virtual void wpgrad(lexer*,fdm*);
	virtual void fillapu(lexer*,fdm*);
	virtual void fillapv(lexer*,fdm*);
	virtual void fillapw(lexer*,fdm*);
	virtual void ptimesave(lexer*,fdm*,ghostcell*);
	void presscorr(lexer*,fdm*,field&);
	void convergence(lexer*,fdm*,ghostcell*);
	void ini_bcval(lexer*,fdm*);

	field1 apu;
    field2 apv;
    field3 apw;
    field4 pcorr;

private:
	int gcval_press, gcval_pcorr;
	int count,q;
	double starttime, endtime;
	int gcval_u, gcval_v, gcval_w;
	double uc,vc,wc,pc;
};

#endif
