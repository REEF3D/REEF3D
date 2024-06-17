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

#ifndef FLUX_FACE_FOU_VRANS_H_
#define FLUX_FACE_FOU_VRANS_H_

#include"flux.h"
#include"increment.h"

class lexer;

using namespace std;

class flux_face_FOU_vrans : public flux, public increment
{
public:

	flux_face_FOU_vrans (lexer *p);
	virtual ~flux_face_FOU_vrans();

	virtual void u_flux(fdm* a,int,field&,double&,double&);
	virtual void v_flux(fdm* a,int,field&,double&,double&);
	virtual void w_flux(fdm* a,int,field&,double&,double&);

private:
    lexer *p;
};

#endif
