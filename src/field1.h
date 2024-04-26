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

#include"field.h"
#include"increment.h"

#ifndef FIELD1_H_
#define FIELD1_H_

using namespace std;

class field1 : public field, public increment
{
public:

	field1 (lexer*);
	virtual ~field1();

    double& operator()(int, int , int);
	double& operator[](int);
    virtual void ggcpol(lexer*);
    virtual void resize(lexer*);
    virtual void dealloc(lexer*);

	int di,dj,dk;
	int imin,imax,jmax,jmin,kmin,kmax;

private:

	void fieldalloc(lexer *);
	void fieldgcalloc(lexer*);
	void fieldlength(lexer *);

	double*** gcfeld;
	int iter;
	int gcfeldsize,feldsize;
	
	int rank, gcextra;
	
	lexer *pp;
};

#endif
