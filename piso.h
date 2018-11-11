/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#ifndef PISO_H_
#define PISO_H_


class piso : public pressure,  public density
{

public:
	piso(lexer* p, fdm *a);
	virtual ~piso();

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
	void presscorr(lexer*,fdm*,field&);
	virtual void ptimesave(lexer*,fdm*,ghostcell*);

	virtual void ucorr2(fdm*,lexer*p,field&);
	virtual void vcorr2(fdm*,lexer*p,field&);
	virtual void wcorr2(fdm*,lexer*p,field&);

	void filluct(lexer*,fdm*);
	void fillvct(lexer*,fdm*);
	void fillwct(lexer*,fdm*);
	
	void ini_bcval(lexer*,fdm*);

	void convergence(lexer*,fdm*,ghostcell*);

    field1 apu,uc;
	field2 apv,vc;
	field3 apw,wc;
	field4 pcorr;

private:
	int gcval_press, gcval_pcorr;
	double starttime, endtime, starttime2, endtime2;
	int cc,q,count,pisoiter;
	double ucc,vcc,wcc,pcc;
	int gcval_u, gcval_v, gcval_w;
};

#endif

