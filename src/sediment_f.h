/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"sediment.h"
#include"sliceint4.h"
#include"slice4.h"
#include"field4a.h"
#include"increment.h"
#include"bedslope.h"
#include"bedload_VR.h"
#include"bedload_einstein.h"
#include"bedload_MPM.h"
#include"bedload_MPM.h"
#include"bedload_EF.h"
#include"bedload_void.h"

class bedload;
class sandslide;
class topo_relax;
class bedshear;
class vrans;
class turbulence;
class sediment_fdm;
class bedshear_reduction;

using namespace std;

#ifndef SEDIMENT_F_H_
#define SEDIMENT_F_H_

class sediment_f : public sediment, public bedslope
{
public:
    sediment_f(lexer*,fdm*,ghostcell*,turbulence*);
	virtual ~sediment_f();
    
    // CFD interface
    virtual void start_cfd(lexer*, fdm*, convection*, ghostcell*, ioflow*, topo*, reinitopo*, suspended*);
    virtual void ini_cfd(lexer*,fdm*,ghostcell*);
    
    virtual void update_cfd(lexer*,fdm*,ghostcell*,ioflow*);
    
    // SFLOW interface
    virtual void start_sflow(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&);
    virtual void ini_sflow(lexer*, fdm2D*, ghostcell*);
    
    void update_sflow(lexer*,fdm2D*,ghostcell*,ioflow*);
    // ---
	
    virtual void relax(lexer*,ghostcell*);
	virtual double bedshear_point(lexer*,fdm*,ghostcell*);
	void sediment_algorithm(lexer*, fdm*, convection*, ghostcell*, ioflow*, topo*, reinitopo*, suspended*);
    
    virtual double qbeval(int,int);
    virtual void qbeget(int,int,double);
    
    void fill_bedk(lexer*,fdm*,ghostcell*);
	void bedlevel(lexer*,fdm*,ghostcell*);
	void topo_zh_update(lexer*,fdm*,ghostcell*,sediment_fdm*);
    void volume_calc(lexer*,fdm*,ghostcell*);
	void filter(lexer*,fdm*,ghostcell*,slice&,int,int);
    
    // print
    virtual void print_3D_bedload(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedload(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu_bedload(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedload(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    
	virtual void print_3D_bedshear(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu_bedshear(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu_bedshear(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_bedshear(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_3D_parameter1(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameter1(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu_parameter1(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameter1(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    
    virtual void print_3D_parameter2(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu_parameter2(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu_parameter2(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu_parameter2(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    

private:
    sediment_fdm *s;
    bedload* pbed;
    sandslide *pslide;
    topo_relax *prelax;
    vrans *pvrans;
    bedshear_reduction *preduce;
	
	bedshear *pbedshear;
    
    double starttime;
    
    int volume_token,sedcalc;
    double volume0;
	
};

#endif
