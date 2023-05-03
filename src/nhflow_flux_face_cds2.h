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

#include"nhflow_flux.h"
#include"increment.h"

#ifndef NHFLOW_FLUX_FACE_CDS2_H_
#define NHFLOW_FLUX_FACE_CDS2_H_

using namespace std;

class nhflow_flux_face_cds2 : public nhflow_flux, public increment
{
public:

	nhflow_flux_face_cds2 (lexer *p);
	virtual ~nhflow_flux_face_cds2();

	virtual void u_flux(fdm_nhf*,int,double*,double&,double&);
	virtual void v_flux(fdm_nhf*,int,double*,double&,double&);
	virtual void w_flux(fdm_nhf*,int,double*,double&,double&);
    virtual void omega_flux(lexer*,fdm_nhf*,int,double*,double*,double*,double&,double&);
    
private:
    lexer *p;

};

#endif
