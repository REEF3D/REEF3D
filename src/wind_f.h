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

#ifndef WIND_F_H_
#define WIND_F_H_

#include"increment.h"
#include"wind.h"

class lexer;
class fdm_nhf;
class ghostcell;
class slice;

using namespace std;

class wind_f final : public wind, public increment
{
public:
    wind_f(lexer*);
	virtual ~wind_f();
    
    void wind_forcing_nhf_x(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*,slice&,slice&) override final;
    void wind_forcing_nhf_y(lexer*, fdm_nhf*, ghostcell*, double*, double*, double*,slice&,slice&) override final;
    
    void wind_forcing_fnpf(lexer*, fdm_fnpf*, ghostcell*,slice&,slice&) override final;

    void wind_forcing_ini(lexer*, ghostcell*) override final;
    
private:
    void wind_forcing_drag_coeff_fnpf(lexer*);
    void wind_forcing_drag_coeff_nhflow(lexer*);
    
    double Cd;
    double cosa,sina;
    double Sx,Sy;
    double xs,xe,ys,ye;
    
    double Uref;

};

#endif
