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

#include"sflow_momentum.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"increment.h"

class sflow_convection;
class sflow_fsf;
class sflow_diffusion;
class sflow_boussinesq;

using namespace std;

#ifndef SFLOW_MOMENTUM_RK4_H_
#define SFLOW_MOMENTUM_RK4_H_

class sflow_momentum_RK4 : public sflow_momentum, public increment
{
public:
	sflow_momentum_RK4(lexer*, fdm2D*, sflow_convection*, sflow_diffusion*, sflow_pressure*, solver2D*, solver2D*, ioflow*, sflow_fsf*, sflow_boussinesq *ppbouss);
	virtual ~sflow_momentum_RK4();
	virtual void start(lexer*, fdm2D*, ghostcell*);

    slice1 Prk1,Prk2,Prk3,Prk;
	slice2 Qrk1,Qrk2,Qrk3,Qrk;
    slice4 wrk1,wrk2,wrk3,wrk;
    slice4 etark1,etark2,etark3,etark,etan;

private:
	void irhs(lexer*,fdm2D*,ghostcell*,slice&,double);
	void jrhs(lexer*,fdm2D*,ghostcell*,slice&,double);
	
	int gcval_u, gcval_v, gcval_w;
	int gcval_urk, gcval_vrk, gcval_wrk,gcval_eta;
	double starttime;

	sflow_convection *pconvec;
	sflow_diffusion *pdiff;
	sflow_pressure *ppress;
	solver2D *psolv;
    solver2D *ppoissonsolv;
	ioflow *pflow;
	sflow_fsf *pfsf;
    sflow_boussinesq *pbouss;
};

#endif
