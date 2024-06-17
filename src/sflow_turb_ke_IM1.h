/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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


#ifndef SFLOW_TURB_KE_IM1_H_
#define SFLOW_TURB_KE_IM1_H_

#include"sflow_turb_io.h"
#include"slice4.h"
#include"sliceint4.h"

class sflow_convection;
class sflow_diffusion;

using namespace std;

class sflow_turb_ke_IM1 : public sflow_turb_io
{

public:
    sflow_turb_ke_IM1(lexer*);
	virtual ~sflow_turb_ke_IM1();
    
	virtual void start(lexer*, fdm2D*, ghostcell*, sflow_convection*, sflow_diffusion*, solver2D*, ioflow*);
	virtual void ktimesave(lexer*, fdm2D*, ghostcell*);
	virtual void etimesave(lexer*, fdm2D*, ghostcell*);
    
private:
    void Pk_update(lexer*, fdm2D*, ghostcell*);
    void ustar_update(lexer*, fdm2D*, ghostcell*);
	void kin_source(lexer*, fdm2D*);
	void eps_source(lexer*, fdm2D*);
    void timesource(lexer*, fdm2D*, slice&);
    void eddyvisc(lexer*, fdm2D*, ghostcell*);
    void clearrhs(lexer*, fdm2D*);
    void wall_law_kin(lexer*, fdm2D*);
    void wall_law_eps(lexer*, fdm2D*);

    slice4 kn, en, Pk, S, ustar, cf;
    sliceint4 wallf;
    double const ce1,ce2,sigk,sige,ceg;
    
    int gcval_kin, gcval_eps;
    int count;
    double starttime;
    
    sflow_convection *pconvec;   
    sflow_diffusion *pdiff; 
    
};

#endif
