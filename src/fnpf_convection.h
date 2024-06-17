/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#ifndef FNPF_CONVECTION_H_
#define FNPF_CONVECTION_H_

class lexer;
class field;
class slice;
class sliceint;
class vec;

using namespace std;

class fnpf_convection
{
public:

    virtual double fx(lexer*, field&, double, double)=0;
	virtual double fy(lexer*, field&, double, double)=0;
	virtual double fz(lexer*, field&, double, double)=0;
    
    virtual double sx(lexer*, slice&, double)=0;
	virtual double sy(lexer*, slice&, double)=0;
    virtual double sz(lexer*, double*)=0;

};

#endif



