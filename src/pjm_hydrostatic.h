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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"pressure.h"
#include"increment.h"

class heat;
class concentration;
class density;

using namespace std;

#ifndef PJM_HYDROSTATIC_H_
#define PJM_HYDROSTATIC_H_

class pjm_hydrostatic : public pressure, public increment
{

public:

	pjm_hydrostatic(lexer* p, fdm *a, heat*&, concentration*&);
	virtual ~pjm_hydrostatic();

	virtual void start(fdm*,lexer*, poisson*, solver*, ghostcell*, ioflow*, field&, field&, field&,double);
    virtual void ini(lexer*,fdm*,ghostcell*);
	virtual void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void vel_setup(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	virtual void ucorr(lexer*,fdm*,field&,double);
	virtual void vcorr(lexer*,fdm*,field&,double);
	virtual void wcorr(lexer*,fdm*,field&,double);
	virtual void upgrad(lexer*,fdm*,slice&,slice&);
	virtual void vpgrad(lexer*,fdm*,slice&,slice&);
    virtual void wpgrad(lexer*,fdm*,slice&,slice&);

private:    
    void debug(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	double starttime,endtime;
	int count, gcval_press;
	int gcval_u, gcval_v, gcval_w;
    
    density *pd;
	
    concentration *pconc;
};



#endif

