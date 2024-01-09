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

class lexer;
class field;
class slice;
class sliceint;
class vec;

using namespace std;

#ifndef FNPF_ETADISC_H_
#define FNPF_ETADISC_H_

class fnpf_etadisc
{
public:

    virtual double sx(lexer*, slice&, slice&)=0;
	virtual double sy(lexer*, slice&, slice&)=0;

};

#endif



