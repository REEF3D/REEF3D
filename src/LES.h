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

#ifndef LES_H_
#define LES_H_

#include"les_io.h"
#include"ghostcell.h"

using namespace std;

class LES : public les_io
{
public:
	LES(lexer *, fdm*);
	virtual ~LES();
	virtual void isource(lexer*, fdm*);
	virtual void jsource(lexer*, fdm*);
	virtual void ksource(lexer*, fdm*);

	double starttime;

};

#endif


