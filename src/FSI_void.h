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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"FSI.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

#ifndef FSI_VOID_H_
#define FSI_VOID_H_

class fsi_void : public fsi
{
public:
	fsi_void(lexer*,ghostcell*){};
	virtual ~fsi_void(){};
	virtual void start(lexer*,fdm*,ghostcell*){};
	virtual void initialize(lexer*,fdm*,ghostcell*){};
    virtual void forcing(lexer*,fdm*,ghostcell*,double,field&,field&,field&,field&,field&,field&,bool){};
    
private:
};

#endif
