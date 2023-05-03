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

#include"printer.h"
#include"increment.h"
#include"field5.h"
class lexer;


#ifndef NODEFILL_H_
#define NODEFILL_H_

using namespace std;

class nodefill : public increment
{

public:
	nodefill(lexer*);
	virtual ~nodefill();
	virtual void nodefill4(lexer*,fdm*,ghostcell*,field&,field&);
    virtual void nodefill4a(lexer*,fdm*,ghostcell*,field&,field&);

private:

    int gcval_phi,gcval_phiext;
    double printtime,printtime2;
};

#endif

