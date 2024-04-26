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

#ifndef field5_H_
#define field5_H_

using namespace std;

class field5 :   public field, public increment
{
public:

	field5 (lexer *);
	virtual ~field5();

    double& operator()(int, int , int);
	double& operator[](int);
    void ggcpol(lexer*);
    virtual void resize(lexer*);
    virtual void dealloc(lexer*);

    //mgc1
	static void makemgc(lexer*);
	static void fillmgc(lexer*);

	const int fip;

private:

	static int a,b,c;
	int di,dj,dk;
	void fieldalloc(lexer *);
	void fieldlength(lexer *);

	double* feld;
	double** gcfeld;
	int iter;
	int gcfeldsize,feldsize;

	int imin,imax,jmax,jmin,kmin,kmax;
	
	lexer *pp;
};

#endif






