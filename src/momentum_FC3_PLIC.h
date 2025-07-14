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
Author: Fabian Knoblauch
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
class heat;
class concentration;
class sixdof;
class fsi;
class VOF_PLIC;
class reini;
class picard;

using namespace std;

#ifndef MOMENTUM_FC3_PLIC
#define MOMENTUM_FC3_PLIC

class momentum_FC3_PLIC : public momentum, public momentum_forcing, public bcmom
{
    
public:
	momentum_FC3_PLIC(lexer*, fdm*, ghostcell*, convection*, diffusion*, pressure*, poisson*, 
                turbulence*, solver*, solver*, ioflow*, heat*&, concentration*&, reini*, fsi*);
	virtual ~momentum_FC3_PLIC();
	virtual void start(lexer*, fdm*, ghostcell*, vrans*,sixdof*,vector<net*>&);
    virtual void utimesave(lexer*, fdm*, ghostcell*);
    virtual void vtimesave(lexer*, fdm*, ghostcell*);
    virtual void wtimesave(lexer*, fdm*, ghostcell*);

    field1 udiff,urk1,urk2,fx;
	field2 vdiff,vrk1,vrk2,fy;
	field3 wdiff,wrk1,wrk2,fz;
    field4 VoF, vof_rk1, vof_rk2;
    
private:
    fluid_update *pupdate;
    
	void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
    void clear_FGH(lexer*,fdm*);
	
	int gcval_u, gcval_v, gcval_w;
    int gcval_vof, gcval_phi,gcval_ro,gcval_visc;
	double starttime;

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
    sixdof *p6dof;
    fsi *pfsi;
    VOF_PLIC *pplic;
    reini *preini;
    density *pd;
    picard *ppicard;
    
};

#endif