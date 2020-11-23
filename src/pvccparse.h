/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef PVCCPARSE_H_
#define PVCCPARSE_H_

class pvccparse : public increment
{
public:
    pvccparse(lexer*,fdm*,ghostcell*);
	virtual ~pvccparse();

	virtual void start(lexer*, fdm*, ghostcell *pgc);
	virtual void ini(lexer*, fdm*, ghostcell *pgc);
	virtual void cellnodes(lexer*, fdm*, ghostcell *pgc);
	virtual void pointcheck(lexer*, fdm*, ghostcell *pgc);
	virtual void collectpoints(lexer*, fdm*, ghostcell *pgc);
	
	virtual void cell_tetra(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_reverse_tetra(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_reverse_tetra_a(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_reverse_tetra_b(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_reverse_tetra_c(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_wedge(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_hexahedron(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide4(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide4a(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide4b(lexer*, fdm*, ghostcell *pgc); 
	virtual void cell_divide4c(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide5a(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide5b(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide5c(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide5d(lexer*, fdm*, ghostcell *pgc);
	virtual void cell_divide6(lexer*, fdm*, ghostcell *pgc);
	
	virtual void face_tetra(lexer*, fdm*, ghostcell *pgc);
	virtual void face_wedge(lexer*, fdm*, ghostcell *pgc);
	virtual void face_hexahedron(lexer*, fdm*, ghostcell *pgc);

private:
    double cx,cy,cz;
    double px,py,pz;
    int pt[16];
    int cl[16];
    int nd[16];
	int fd[6][20],fcount[6];
    const double dx,eps;
    int pcount,nn,cc,clcount,pointcount,count,n,q;
	int pcount0,pcount1;
    int check1,check2,check3,check4,check9,check10,check11,check12;
};

#endif

