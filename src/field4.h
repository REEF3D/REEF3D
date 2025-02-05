/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef FIELD4_H_
#define FIELD4_H_

#include"field.h"
#include"increment.h"

using namespace std;

class field4 : public field, increment
{
public:

	field4 (lexer*);
	virtual ~field4();

    double& operator()(int, int , int);
	double& operator[](int);
    virtual void resize(lexer*);
    virtual void dealloc(lexer*);
    
	int di,dj,dk;
	int imin,imax,jmax,jmin,kmin,kmax;

private:

	void fieldalloc(lexer *);
	void fieldlength(lexer *);

    int iter;
	int feldsize;
    
    int level[12];
	
	int rank, gcextra;
	
	double starttime;
	
	lexer *pp;

};

#endif





