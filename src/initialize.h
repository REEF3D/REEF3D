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

#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include"resize.h"
#include"increment.h"

class fdm;
class lexer;
class ghostcell;
class turbulence;
class sediment;

using namespace std;

class initialize : public increment, private resize_class
{

public:
	initialize(lexer*);
	virtual ~initialize();

	void start(fdm*, lexer*, ghostcell*);
    void droplet_ini(lexer*,fdm*,ghostcell*);
	void hydrostatic(lexer*,fdm*,ghostcell*);
	void iniphi_io(fdm*, lexer*,ghostcell*);
	void inivof_io(fdm*, lexer*,ghostcell*);
	void iniphi_surfarea(lexer*,fdm*,ghostcell*);
	void stateini(lexer*,fdm*,ghostcell*,turbulence*,sediment*);
    void inipsi(lexer*,fdm*,ghostcell*);

private:
	void inifdm(fdm*, lexer*, ghostcell*);
	void iniphi(fdm*, lexer*,ghostcell*);
	void iniphi_box(lexer*,fdm*,ghostcell*);	
	void inivof(fdm*, lexer*,ghostcell*);
	void inivof_box(lexer*,fdm*,ghostcell*);
	void inivofPLIC(fdm*, lexer*,ghostcell*);
	void bcwall_check(fdm*,lexer*);
	void nodecalc(fdm*, lexer*);
	void faceneighbors(lexer*,fdm*);
	void maxcoor(fdm*, lexer*,ghostcell*);
	void paraini(lexer*, fdm*,ghostcell*);
	void pressini(lexer*,fdm*,ghostcell*);
	void topoini(lexer*,fdm*,ghostcell*);
    
	int conv(double);

	const double smallnum;
	double epsi;

	int n,q,iend,kend;
	double deltax;
	double H;
};

#endif
