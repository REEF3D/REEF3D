/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

class fdm;
class lexer;
class field;

#ifndef FILLVEC_H_
#define FILLVEC_H_

using namespace std;

class fillvec : virtual public increment
{
public:

	fillvec ();
	virtual ~fillvec();
	
	virtual void fillxvec1(lexer*,fdm*,field&);
    virtual void fillxvec2(lexer*,fdm*,field&);
    virtual void fillxvec3(lexer*,fdm*,field&);
    virtual void fillxvec4(lexer*,fdm*,field&);

	virtual void fillxfield1(lexer*,fdm*,field&);
	virtual void fillxfield2(lexer*,fdm*,field&);	
	virtual void fillxfield3(lexer*,fdm*,field&);
	virtual void fillxfield4(lexer*,fdm*,field&);
	

private:
	int margin;
	int count,q;
	
};

#endif
