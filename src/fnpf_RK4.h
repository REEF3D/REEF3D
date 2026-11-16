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

#ifndef FNPF_RK4_H_
#define FNPF_RK4_H_

#include"fnpf_rk.h"

using namespace std;

// FNPF time step with the classical fourth-order Runge-Kutta scheme, A 310 4
class fnpf_RK4 final : public fnpf_rk
{
public:
	fnpf_RK4(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_RK4();

    void start(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*) override final;
    void ini_wetdry(lexer*, fdm_fnpf*, ghostcell*) override final;

private:
    slice4 erk1,erk2,erk3,erk4;   // stage tendencies deta/dt
    slice4 frk1,frk2,frk3,frk4;   // stage tendencies dFifsf/dt
    slice4 erk,frk;               // stage values
};

#endif
