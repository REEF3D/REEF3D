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

#include"increment.h"
#include"fieldint.h"

class lexer;

#ifndef INTFIELD_H_
#define INTFIELD_H_

using namespace std;

class fieldint5 :  public fieldint, public increment
{
public:

	fieldint5 (lexer *);
	virtual ~fieldint5();

    int& operator()(int, int , int);
    
    virtual void resize(lexer*);

private:

	int di,dj,dk;
	void fieldalloc(lexer *);
	void fieldlength(lexer *);
	void fieldgcalloc(lexer*);

	int* feld;
	int iter;

	static int imin,imax,jmax,jmin,kmin,kmax;
	
	lexer *pp;

};

#endif






