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

#ifndef NHFLOW_SCALAR_ADVEC_H_
#define NHFLOW_SCALAR_ADVEC_H_

class lexer;
class fdm_nhf;

using namespace std;

class nhflow_scalar_advec
{
public:

    virtual void uadvec(int,double*,double&,double&)=0;
	virtual void vadvec(int,double*,double&,double&)=0;
	virtual void wadvec(int,double*,double&,double&)=0;

};

#endif
