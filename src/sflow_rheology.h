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

class lexer;
class fdm2D;
class slice;

#ifndef SFLOW_RHEOLOGY_H_
#define SFLOW_RHEOLOGY_H_

using namespace std;

class sflow_rheology
{

public:

	virtual void u_source(lexer*, fdm2D*, slice&, slice&)=0;
    virtual void v_source(lexer*, fdm2D*, slice&, slice&)=0;

};

#endif
