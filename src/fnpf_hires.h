/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#ifndef FNPF_HIRES_H_
#define FNPF_HIRES_H_

#include"fnpf_convection.h"
#include"increment.h"

class fdm_fnpf;

using namespace std;

class fnpf_hires : public fnpf_convection, public increment
{
public:
    fnpf_hires(lexer*,fdm_fnpf*);
	virtual ~fnpf_hires();

    virtual double fx(lexer*, field&, double, double){};
	virtual double fy(lexer*, field&, double, double){};
	virtual double fz(lexer*, field&, double, double){};
    
    virtual double sx(lexer*, slice&, double);
	virtual double sy(lexer*, slice&, double);
    virtual double sz(lexer*, double*){};
    
private:
    fdm_fnpf *c;
    
    double limiter(double v1, double v2);
    
    double dfdx_min, dfdx_plus, dfdy_min, dfdy_plus, dfdz_min, dfdz_plus;
    double denom,val,grad;

};

#endif
