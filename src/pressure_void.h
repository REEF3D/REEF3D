/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#ifndef PRESSURE_VOID_H_
#define PRESSURE_VOID_H_

#include"pressure.h"
#include"increment.h"

using namespace std;

class pressure_void final : public pressure, public increment
{

public:
	pressure_void(lexer* p);
	virtual ~pressure_void();

	void start(fdm*,lexer* p, poisson*, solver*, ghostcell*,ioflow*, field&, field&, field&,double) override final;
    void ini(lexer*,fdm*,ghostcell*) override final;
	void rhs(lexer*,fdm*,ghostcell*,field&,field&,field&,double);
	void ucorr(lexer*p,fdm*,field&,double) override final;
	void vcorr(lexer*p,fdm*,field&,double) override final;
	void wcorr(lexer*p,fdm*,field&,double) override final;
	void upgrad(lexer*,fdm*,slice&,slice&) override final;
	void vpgrad(lexer*,fdm*,slice&,slice&) override final;
	void wpgrad(lexer*,fdm*,slice&,slice&) override final;
};

#endif
