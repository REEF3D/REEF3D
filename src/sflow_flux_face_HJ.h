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

#ifndef SFLOW_FLUX_FACE_HJ_H_
#define SFLOW_FLUX_FACE_HJ_H_

#include"increment.h"
#include"sflow_flux.h"

class lexer;
class fdm2D;

using namespace std;

class sflow_flux_face_HJ : public sflow_flux, public increment
{
public:

	sflow_flux_face_HJ (lexer *p);
	virtual ~sflow_flux_face_HJ();

	virtual void u_flux(int,slice&,double&,double&);
	virtual void v_flux(int,slice&,double&,double&);

private:
    lexer *p;
    
};

#endif
