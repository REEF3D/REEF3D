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

#include"sediment.h"
#include"topo_vel.h"
#include"sliceint4.h"
#include"slice4.h"
#include"field4a.h"

class sandslide;
class topo_relax;
class bedshear;

using namespace std;

#ifndef SEDIMENT_F_H_
#define SEDIMENT_F_H_

class sediment_f : public sediment, topo_vel
{
public:
    sediment_f(lexer*,turbulence*);
	virtual ~sediment_f();

	virtual void start(lexer*, fdm*, convection*, ghostcell*, ioflow*, topo*, reinitopo*, suspended*, bedload*);
	virtual void update(lexer*,fdm*,ghostcell*,ioflow*);
    virtual void relax(lexer*,fdm*,ghostcell*);
    virtual void ini(lexer*,fdm*,ghostcell*);
	virtual double bedshear_point(lexer*,fdm*,ghostcell*);
	void sediment_algorithm(lexer*, fdm*, convection*, ghostcell*, ioflow*, topo*, reinitopo*, suspended*, bedload*);
    
	void fill_bss(lexer*,fdm*,ghostcell*);
    void fill_bedk(lexer*,fdm*,ghostcell*);
	void timestep(lexer*,fdm*,ghostcell*);
	void bedlevel(lexer*,fdm*,ghostcell*);
	void topo_zh_update(lexer*,fdm*,ghostcell*);
	void filter(lexer*,fdm*,ghostcell*,slice&,int,int);
    
	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    

private:
    sandslide *pslide;
    topo_relax *prelax;
	
	bedshear *pbedshear;
    
    field4a zh,bss;
    slice4 bedtau;
    
    double starttime;
	
};

#endif
