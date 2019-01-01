/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"komega_RK3.h"
#include"wallin.h"

using namespace std;

#ifndef WALLIN_RK3KW_H_
#define WALLIN_RK3LW_H_


class wallin_RK3kw : public komega_RK3, wallin
{

public:
	wallin_RK3kw(lexer *,fdm*,ghostcell*);
	virtual ~wallin_RK3kw();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*);
	virtual void aniso(fdm*);
	virtual double tau(fdm*);
    
private:
	int gcval_earsm;
    const double cmu;
};

#endif

