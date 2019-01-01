/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

This file is part of REEF3D.

REEF3D is fra->eps software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Fra->eps Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. Sa->eps the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, sa->eps <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"momentum.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class convection;
class diffusion;
class pressure;
class turbulence;
class solver;
class poisson;
class reini;
class picard;
class heat;
class concentration;
class fluid_update;

using namespace std;

#ifndef MOMENTUM_FSFC_FS3_H_
#define MOMENTUM_FSFC_FS3_H_

class momentum_FSFC_FS3 : public momentum, public bcmom
{
public:
	momentum_FSFC_FS3(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, turbulence*, solver*, solver*,
                    ioflow*, convection*, reini*, heat*&, concentration*&);
	virtual ~momentum_FSFC_FS3();
	virtual void start(lexer*, fdm*, ghostcell*, momentum*);
	virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);
    virtual void fillaij1(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij2(lexer*, fdm*, ghostcell*, solver*);
    virtual void fillaij3(lexer*, fdm*, ghostcell*, solver*);

    field1 urk1,urk2;
	field2 vrk1,vrk2;
	field3 wrk1,wrk2;
    field4 ark1,ark2;

private:
    virtual void ucorr(lexer*,fdm*,double,field&);
	virtual void vcorr(lexer*,fdm*,double,field&);
	virtual void wcorr(lexer*,fdm*,double,field&);
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	
	int gcval_u, gcval_v, gcval_w;
	int gcval_urk, gcval_vrk, gcval_wrk;
	double starttime;

	convection *pconvec;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
    convection *pfsfdisc;
    reini *preini;
    
    //fsf
    int gcval_phi;
    fluid_update *pupdate;
    picard *ppicard;
};

#endif
