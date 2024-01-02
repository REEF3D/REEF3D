/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nsewave.h"
#include"bcmom.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"

class picard;
class fluid_update;
class heat;
class concentration;

using namespace std;

#ifndef NSEWAVE_RK3_H_
#define NSEWAVE_RK3_H_

class nsewave_RK3 : public nsewave, public bcmom
{
public:
    nsewave_RK3(lexer*, fdm*, ghostcell*,heat*&,concentration*&);
	virtual ~nsewave_RK3();
    
    virtual void start(lexer*, fdm*, ghostcell*, momentum*, diffusion*, turbulence*, convection*, pressure*, 
                        poisson*, solver*, solver*, ioflow*, vrans*, sixdof_df_base*, vector<net*>&);
    virtual void ini(lexer*, fdm*, ghostcell*, ioflow*);
    void update(lexer*,fdm*,ghostcell*,slice&);
    
private: 

    void eta_disc(lexer*,fdm*,ghostcell*,field&,field&);
    void irhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void jrhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
	void krhs(lexer*,fdm*,ghostcell*,field&,field&,field&,field&,double);
    
    field1 urk1,urk2;
	field2 vrk1,vrk2;
	field3 wrk1,wrk2;
    
    fluid_update *pupdate;
    picard *ppicard;
    
    int gcval_phi;
	double starttime;
    double phival,H;
    const double epsi;
    int gcval_u, gcval_v, gcval_w;
    
    slice1 P;
    slice2 Q;
    slice4 eta,etark1,etark2,L;

};

#endif
