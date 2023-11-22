/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"resize.h"
#include"increment.h"

class fdm_ptf;
class lexer;
class ghostcell;

#ifndef INITIALISE_PTF_H_
#define INITIALISE_PTF_H_

using namespace std;

class initialise_ptf : public increment, private resize_class
{

public:
	initialise_ptf(lexer*);
	virtual ~initialise_ptf();

	void start(fdm_ptf*, lexer*, ghostcell*);
	void hydrostatic(lexer*,fdm_ptf*,ghostcell*);
	void iniphi_io(fdm_ptf*, lexer*,ghostcell*);
	void iniphi_surfarea(lexer*,fdm_ptf*,ghostcell*);
	void stateini(lexer*,fdm_ptf*,ghostcell*);
    void inipsi(lexer*,fdm_ptf*,ghostcell*);

private:
	void inifdm_ptf(fdm_ptf*, lexer*, ghostcell*);
	void iniphi(fdm_ptf*, lexer*,ghostcell*);
	void iniphi_box(lexer*,fdm_ptf*,ghostcell*);	
	void bcwall_check(fdm_ptf*,lexer*);
	void nodecalc(fdm_ptf*, lexer*);
	void faceneighbors(lexer*,fdm_ptf*);
	void maxcoor(fdm_ptf*, lexer*,ghostcell*);
	void paraini(lexer*, fdm_ptf*,ghostcell*);
	void pressini(lexer*,fdm_ptf*,ghostcell*);
	void topoini(lexer*,fdm_ptf*,ghostcell*);
    
	int conv(double);

	const double smallnum;
	double epsi;

	int n,q,iend,kend;
	double deltax;
	double H;
};

#endif
