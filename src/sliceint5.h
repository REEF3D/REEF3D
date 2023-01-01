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

#include"sliceint.h"
#include"increment.h"

#ifndef SLICEINT5_H_
#define SLICEINT5_H_

using namespace std;

class sliceint5 : public sliceint, increment
{
public:

	sliceint5 (lexer*);
	virtual ~sliceint5();

    virtual int& operator()(int, int);

    virtual void resize(lexer*);
    virtual void dealloc(lexer*);
    
	int di,dj;
	int imin,imax,jmax,jmin;

	int ***gcfeld;

private:

	void fieldalloc(lexer *);

    int iter;
	int feldsize;
	
	
	double starttime;
	
	lexer *pp;

};

#endif





