/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

class lexer;
class slice;

#ifndef NHFLOW_HXY_FLUX_H_
#define NHFLOW_HXY_FLUX_H_

using namespace std;

class nhflow_hxy_flux : public increment
{
public:

	nhflow_hxy_flux (lexer *p);
	virtual ~nhflow_hxy_flux();

	void u_flux(int,slice&,double&,double&);
	void v_flux(int,slice&,double&,double&);

private:
    lexer *p;
    
};

#endif
