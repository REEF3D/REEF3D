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
class ghostcell;

#ifndef MGCCART_H_
#define MGCCART_H_

using namespace std;

class mgccart :  public increment
{
public:

	mgccart (lexer *);
	virtual ~mgccart();

    //mgc
	void startmgc(lexer*);
	void makemgc(lexer*,int*);
	void mgcsetup(lexer*,int*);
	void gcdirfill(lexer*,int,int**,int*);
	void fillmgc(lexer*,int,int**,int*);
	
	void gcsidefill(lexer*);

	
	int imin,imax,jmax,jmin,kmin,kmax;
	int gcdirsize;
	
	
};

#endif






