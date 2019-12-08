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

#include"sflow_sediment.h"
#include"increment.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class fnpf_convection;

using namespace std;

#ifndef SFLOW_SEDIMENT_F_H_
#define SFLOW_SEDIMENT_F_H_

class sflow_sediment_f : public sflow_sediment, public increment
{
public:

    sflow_sediment_f(lexer*, fdm2D*);
	virtual ~sflow_sediment_f();

	virtual void ini(lexer*, fdm2D*, ghostcell*);
    virtual void start(lexer*, fdm2D*, ghostcell*, slice&, slice&,slice&);
    
private:
    void sediment_algorithm(lexer*, fdm2D*, ghostcell*,slice&,slice&,slice&);
    
    void bedslope(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void bednode(lexer*, fdm2D*, ghostcell*);
    void bedshear(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void shields(lexer*, fdm2D*, ghostcell*);
    void bedshear_slope(lexer*, fdm2D*, ghostcell*);
    void dey_ana(lexer*,fdm2D*,ghostcell*);
    void dey_emp(lexer*,fdm2D*,ghostcell*);
    void parker(lexer*,fdm2D*,ghostcell*);
    void fredsoe_long(lexer*,fdm2D*,ghostcell*);
    
    void bedload(lexer*, fdm2D*, ghostcell*);
    void bedload_vanRijn(lexer*, fdm2D*, ghostcell*);
    
    void exner(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&);
    void bedchange_update(lexer*, fdm2D*, ghostcell*);
    
    void sandslide(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void sandslide_v2(lexer*, fdm2D*, ghostcell*, slice&, slice&);
    void slide(lexer*, fdm2D*, ghostcell*);
    void slide_v2(lexer*, fdm2D*, ghostcell*);
    
    void filter(lexer*, fdm2D*, ghostcell*,slice&,int,int);
    
    
    slice4 tau,taucr,alpha,teta,gamma,phi,fh,red;
    slice4 topovel1,topovel2;
    
    double starttime;
    
    fnpf_convection *pdx;
    
    double delta, midphi;
    double fac1,fac2;
    int slidecount;
    
    // relax
    void relax_ini(lexer*, fdm2D*);
    void relax(lexer*, fdm2D*,ghostcell*);
    virtual double rf(lexer*, fdm2D*,ghostcell*);

	double distcalc(lexer*,double, double, double);
	double r1(lexer*, double, double);
	
	double *tan_betaS73,*betaS73,*dist_S73;
	double val;
    
};

#endif