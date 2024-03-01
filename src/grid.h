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

#include"increment.h"

class lexer;
class fdm;
class ghostcell;
class fieldint;
class cpt;

#ifndef GRID_H_
#define GRID_H_

using namespace std;

class grid :  public increment
{
public:

	grid (lexer *);
	virtual ~grid();
	
	void makegrid(lexer*,ghostcell*);
	void update_topo_grid(lexer*,ghostcell*);
	void update_sixdof_grid(lexer*,ghostcell*);
    
    // CPT
    void column_pt1_update(lexer*,cpt&);
    void column_pt2_update(lexer*,cpt&);
    void column_pt3_update(lexer*,cpt&);
    void column_pt4_update(lexer*,cpt&);
    void column_pt4a_update(lexer*,cpt&);
    void column_pt6_update(lexer*,cpt&);

	void column_pt1_assign(lexer*,fieldint&,cpt&);
	void column_pt2_assign(lexer*,fieldint&,cpt&);
	void column_pt3_assign(lexer*,fieldint&,cpt&);
	void column_pt4_assign(lexer*,fieldint&,cpt&);
    void column_pt4a_assign(lexer*,fieldint&,cpt&);
    void column_pt6_assign(lexer*,fieldint&,cpt&);

    int column_pt1_count(lexer*);
	int column_pt2_count(lexer*);
	int column_pt3_count(lexer*);
	int column_pt4_count(lexer*);
    int column_pt4a_count(lexer*);
    int column_pt6_count(lexer*);
    
    
    // cval
    void cval_update1(lexer*,fieldint&);
	void cval_update2(lexer*,fieldint&);
	void cval_update3(lexer*,fieldint&);
	void cval_update4(lexer*,fieldint&);
    void cval_update4a(lexer*,fieldint&);
    void cval_update6(lexer*,fieldint&);

	void cval_gcb1(lexer*,fieldint&);
    void cval_gcb2(lexer*,fieldint&);
    void cval_gcb3(lexer*,fieldint&);
    void cval_gcb4(lexer*,fieldint&);
    void cval_gcb4a(lexer*,fieldint&);
    void cval_gcb6(lexer*,fieldint&);

	void cval_gcpara1(lexer*,fieldint&);
    void cval_gcpara2(lexer*,fieldint&);
    void cval_gcpara3(lexer*,fieldint&);
    void cval_gcpara4(lexer*,fieldint&);
    void cval_gcpara4a(lexer*,fieldint&);
    void cval_gcpara6(lexer*,fieldint&);

    
private:

    int g,q,margin,count;

   
	
};

#endif






