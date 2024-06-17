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

#ifndef NHFLOW_FLUX_BUILD_F_H_
#define NHFLOW_FLUX_BUILD_F_H_

#include"nhflow_flux_build.h"
#include"increment.h"

class patchBC_interface;

using namespace std;

class nhflow_flux_build_f : public nhflow_flux_build, public increment
{

public:

	nhflow_flux_build_f(lexer*,ghostcell*,patchBC_interface*);
	virtual ~nhflow_flux_build_f();

    virtual void start_E(lexer*, fdm_nhf*, ghostcell*);
    virtual void start_U(lexer*, fdm_nhf*, ghostcell*);
    virtual void start_V(lexer*, fdm_nhf*, ghostcell*);
    virtual void start_W(lexer*, fdm_nhf*, ghostcell*);

private:
    patchBC_interface *pBC;

};

#endif
