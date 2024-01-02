/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MEFCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"momentum.h"
#include"momentum_forcing.h"
#include"bcmom.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class convection;
class diffusion;
class pressure;
class density;
class turbulence;
class solver;
class poisson;
class fluid_update;
class nhflow;
class reini;
class picard;
class heat;
class concentration;
class sixdof_df_base;
class fsi;

using namespace std;

#ifndef MOMENTUM_FCC3_H_
#define MOMENTUM_FCC3_H_

class momentum_FCC3 : public momentum, public momentum_forcing, public bcmom
{
public:
	momentum_FCC3(lexer*, fdm*, ghostcell*, convection*, convection*, diffusion*, pressure*, poisson*, 
                turbulence*, solver*, solver*, ioflow*, heat*&, concentration*&, reini*, fsi*);
	virtual ~momentum_FCC3();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof_df_base*,vector<net*>&);
    virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);

    field1 ur,udiff,urk1,urk2,fx;
	field2 vr,vdiff,vrk1,vrk2,fy;
	field3 wr,wdiff,wrk1,wrk2,fz;
    
    field1 Mx,rox;
    field2 My,roy;
    field3 Mz,roz;
    
    field1 Mx_rk1,Mx_rk2;
	field2 My_rk1,My_rk2;
	field3 Mz_rk1,Mz_rk2;
    
    field1 rox_rk1,rox_rk2;
	field2 roy_rk1,roy_rk2;
	field3 roz_rk1,roz_rk2;
    
    field4 ls,frk1,frk2;

private:
    fluid_update *pupdate;
    picard *ppicard;
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	
    void clear_FGH(lexer*,fdm*);
    void face_density(lexer*,fdm*,ghostcell*,field&,field&,field&);
    
    
    double vel_limiter(lexer*,fdm*,field&,field&,field&,field&);
    double ro_filter(lexer*,fdm*,field&);

    
	int gcval_u, gcval_v, gcval_w;
    int gcval_phi;
    double val;
	double starttime;
    double ro_threshold;

	convection *pconvec;
    convection *pfsfdisc;
	diffusion *pdiff;
	pressure *ppress;
	poisson *ppois;
	turbulence *pturb;
	solver *psolv;
    solver *ppoissonsolv;
	ioflow *pflow;
    nhflow *pnh;
    reini *preini;
    density *pd;
    sixdof_df_base *p6dof_df;
    fsi *pfsi;
    
};

#endif
