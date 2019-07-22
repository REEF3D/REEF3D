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
#include"increment.h"

class density;
class heat;
class concentration;

using namespace std;

#ifndef PJM_SIG_H_
#define PJM_SIG_H_

class pjm_sig : public pressure, public increment
{

public:

	pjm_sig(lexer* p, fdm *a, heat*&, concentration*&);
	virtual ~pjm_sig();

	virtual void start(fdm*,lexer* p, poisson*, solver*, ghostcell*,momentum*,ioflow*, field&, field&, field&,double);
	virtual void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void vel_setup(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void ucorr(lexer*p,fdm*,field&,double);
	virtual void vcorr(lexer*p,fdm*,field&,double);
	virtual void wcorr(lexer*p,fdm*,field&,double);
	virtual void upgrad(lexer*,fdm*);
	virtual void vpgrad(lexer*,fdm*);
	virtual void ptimesave(lexer*,fdm*,ghostcell*);

	virtual void wpgrad(lexer*,fdm*);
	virtual void fillapu(lexer*,fdm*);
	virtual void fillapv(lexer*,fdm*);
	virtual void fillapw(lexer*,fdm*);
	virtual void pressure_norm(lexer*,fdm*,ghostcell*);


private:
	double starttime,endtime;
	int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
	
	void debug(lexer*,fdm*);
    
    density *pd;
};



#endif

