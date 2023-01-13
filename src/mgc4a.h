/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#ifndef MGC4A_H_
#define MGC4A_H_

using namespace std;

class mgc4a :  public increment
{
public:

	mgc4a (lexer *);
	virtual ~mgc4a();

    //mgc4
	void makemgc(lexer*);
	void mgcsetup(lexer*);
	void gcdirfill(lexer*);
	void fillmgc(lexer*);
    void fillgcb(lexer*);

	// ggc
	void make_ggc(lexer*);
	void fill_ggc(lexer*);
    
    // dgc
    void make_dgc(lexer*);
    void fill_dgc(lexer*);

	int imin,imax,jmax,jmin,kmin,kmax;
    
private:
	int di,dj,dk;
	int qn;
    
    int *hgc;
};

#endif







