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

#ifndef SLICE_H_
#define SLICE_H_

class lexer;
class fdm;

using namespace std;

class slice
{
public:
	virtual double& operator()(int, int)=0;
	virtual double& operator[](int)=0;
	virtual void ggcpol(lexer*)=0;
    virtual void resize(lexer*)=0;
    virtual void dealloc(lexer*)=0;
	
	double *V;
};

#endif






