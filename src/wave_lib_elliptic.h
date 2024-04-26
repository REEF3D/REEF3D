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

class lexer;
class fdm;
class ghostcell;

#ifndef WAVE_LIB_ELLIPTIC_H_
#define WAVE_LIB_ELLIPTIC_H_

using namespace std;

class wave_lib_elliptic 
{
public:
    wave_lib_elliptic();
	virtual ~wave_lib_elliptic();
    
    void elliptic(lexer*,double,double&,double&,double&);
	double K_elliptic_1(double);
	double E_elliptic_1(double);
	double K_elliptic_5(double);
	double E_elliptic_5(double);
	double K_elliptic(double);
	double E_elliptic(double);
    
    double Km,Em,ell,eta2;
    const double epsi;
    double modulus;
    

};

#endif
