/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"komega_IM2.h"
#include"EARSM.h"

class multiphase;

using namespace std;

#ifndef EARSM_KW_IM2_H_
#define EARSM_KW_IM2_H_


class EARSM_kw_IM2 : public komega_IM2, public EARSM
{

public:
	EARSM_kw_IM2(lexer *,fdm*,ghostcell*);
	virtual ~EARSM_kw_IM2();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*);
	virtual void aniso(fdm*);
	virtual double tau(fdm*);
    
private:
	int gcval_earsm;
    const double cmu;
};

#endif

