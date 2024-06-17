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

#ifndef MGCSLICE4_H_
#define MGCSLICE4_H_

#include"increment.h"

class lexer;
class ghostcell;

using namespace std;

class mgcslice4 :  public increment
{
public:

	mgcslice4 (lexer *);
	virtual ~mgcslice4();

    //mgcslice4
	void makemgc(lexer*);
	void mgcsetup(lexer*);
	void gcdirfill(lexer*);
	void fillmgc(lexer*);
	void make_ggc(lexer*);
	void fill_ggc(lexer*);

    
    void gcb_seed(lexer*);
	
	int imin,imax,jmax,jmin,kmin,kmax;
	int gcdirsize;
	int ggcsize;
	
private:
	int di,dj;
	int qn;
    int count;
	
};

#endif






